#include "ovdbutil/hollowing.h"
#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetRebuild.h>

#include "trimesh2/TriMesh.h"
#include "mmesh/trimesh/trimeshutil.h"

#include "ccglobal/tracer.h"
#include <openvdb/tools/RayIntersector.h>

namespace ovdbutil
{
    inline trimesh::vec3 to_vec3f(const openvdb::Vec3s& v) { return trimesh::vec3(v.x(), v.y(), v.z()); }
    inline trimesh::dvec3 to_vec3d(const openvdb::Vec3s& v) { return trimesh::dvec3(v.x(), v.y(), v.z()); }
    inline trimesh::ivec3 to_vec3i(const openvdb::Vec3I& v) { return trimesh::ivec3(int(v[0]), int(v[1]), int(v[2])); }
    inline trimesh::ivec4 to_vec4i(const openvdb::Vec4I& v) { return trimesh::ivec4(int(v[0]), int(v[1]), int(v[2]), int(v[3])); }

    struct Contour3D {
        std::vector<trimesh::vec3> points;
        std::vector<trimesh::ivec3> faces3;
        std::vector<trimesh::ivec4> faces4;

        Contour3D() = default;
        inline bool empty() const
        {
            return points.empty() || (faces4.empty() && faces3.empty());
        }
    };

    trimesh::TriMesh* to_triangle_mesh(const Contour3D& ctour) {
        trimesh::TriMesh* mesh = new trimesh::TriMesh();
        mesh->vertices = ctour.points;
        if (ctour.faces4.empty())
        {
            mesh->faces = ctour.faces3;
        }
        else
        {
            mesh->faces.reserve(ctour.faces3.size() + 2 * ctour.faces4.size());
            std::copy(ctour.faces3.begin(), ctour.faces3.end(),
                std::back_inserter(mesh->faces));

            for (const trimesh::ivec4& quad : ctour.faces4) {
                mesh->faces.emplace_back(quad[0], quad[1], quad[2]);
                mesh->faces.emplace_back(quad[2], quad[3], quad[0]);
            }
        }

        return mesh;
    }

    class TriangleMeshDataAdapter {
    public:
        trimesh::TriMesh& mesh;

        size_t polygonCount() const { return mesh.faces.size(); }
        size_t pointCount() const { return mesh.vertices.size(); }
        size_t vertexCount(size_t) const { return 3; }

        // Return position pos in local grid index space for polygon n and vertex v
        void getIndexSpacePoint(size_t n, size_t v, openvdb::Vec3d& pos) const;
    };

    class Contour3DDataAdapter {
    public:
        const Contour3D& mesh;

        size_t polygonCount() const { return mesh.faces3.size() + mesh.faces4.size(); }
        size_t pointCount() const { return mesh.points.size(); }
        size_t vertexCount(size_t n) const { return n < mesh.faces3.size() ? 3 : 4; }

        // Return position pos in local grid index space for polygon n and vertex v
        void getIndexSpacePoint(size_t n, size_t v, openvdb::Vec3d& pos) const;
    };

    void TriangleMeshDataAdapter::getIndexSpacePoint(size_t          n,
        size_t          v,
        openvdb::Vec3d& pos) const
    {
        size_t vidx = size_t(mesh.faces[n][v]);
        trimesh::dvec3 p = trimesh::dvec3(mesh.vertices[vidx]);
        pos = { p.x, p.y, p.z };
    }

    void Contour3DDataAdapter::getIndexSpacePoint(size_t          n,
        size_t          v,
        openvdb::Vec3d& pos) const
    {
        size_t vidx = 0;
        if (n < mesh.faces3.size()) vidx = size_t(mesh.faces3[n][v]);
        else vidx = size_t(mesh.faces4[n - mesh.faces3.size()][v]);

        trimesh::vec3 p = mesh.points[vidx];
        pos = { p.x, p.y, p.z };
    }

// TODO: Do I need to call initialize? Seems to work without it as well but the
// docs say it should be called ones. It does a mutex lock-unlock sequence all
// even if was called previously.
    openvdb::FloatGrid::Ptr mesh_to_grid(trimesh::TriMesh* mesh,
        const openvdb::math::Transform& tr,
        float               exteriorBandWidth,
        float               interiorBandWidth,
        double voxel_size,
        int                 flags = 0
        )
    {
        openvdb::initialize();

        openvdb::FloatGrid::Ptr grid = openvdb::tools::meshToVolume<openvdb::FloatGrid>(
                TriangleMeshDataAdapter{ *mesh }, tr, exteriorBandWidth,
                interiorBandWidth, voxel_size, flags);

        return grid;
    }

    openvdb::FloatGrid::Ptr mesh_to_grid(const Contour3D& mesh,
        const openvdb::math::Transform& tr,
        float exteriorBandWidth,
        float interiorBandWidth,
        int flags = 0)
    {
        openvdb::initialize();
        return openvdb::tools::meshToVolume<openvdb::FloatGrid>(
            Contour3DDataAdapter{ mesh }, tr, exteriorBandWidth, interiorBandWidth,
            flags);
    }

    template<class Grid>
    Contour3D _volumeToMesh(const Grid& grid,
        double      isovalue,
        double      adaptivity,
        bool        relaxDisorientedTriangles)
    {
        openvdb::initialize();

        std::vector<openvdb::Vec3s> points;
        std::vector<openvdb::Vec3I> triangles;
        std::vector<openvdb::Vec4I> quads;

        openvdb::tools::volumeToMesh(grid, points, triangles, quads, isovalue,
            adaptivity, relaxDisorientedTriangles);

        Contour3D ret;
        ret.points.reserve(points.size());
        ret.faces3.reserve(triangles.size());
        ret.faces4.reserve(quads.size());

        for (auto& v : points) ret.points.emplace_back(to_vec3d(v));
        for (auto& v : triangles) ret.faces3.emplace_back(to_vec3i(v));
        for (auto& v : quads) ret.faces4.emplace_back(to_vec4i(v));

        return ret;
    }

    trimesh::TriMesh* grid_to_mesh(const openvdb::FloatGrid& grid,
        double                    isovalue,
        double                    adaptivity,
        bool                      relaxDisorientedTriangles)
    {
        return to_triangle_mesh(
            _volumeToMesh(grid, isovalue, adaptivity, relaxDisorientedTriangles));
    }

    Contour3D grid_to_contour3d(const openvdb::FloatGrid& grid,
        double                    isovalue,
        double                    adaptivity,
        bool relaxDisorientedTriangles)
    {
        return _volumeToMesh(grid, isovalue, adaptivity,
            relaxDisorientedTriangles);
    }

    openvdb::FloatGrid::Ptr redistance_grid(const openvdb::FloatGrid& grid, double iso, double er = 3.0, double ir = 3.0)
    {
        return openvdb::tools::levelSetRebuild(grid, float(iso), float(er), float(ir));
    }

	double Cross(const trimesh::vec2& p1, const trimesh::vec2& p2, const trimesh::vec2& p3, const trimesh::vec2& p4)
	{
		return (p2.x - p1.x) * (p4.y - p3.y) - (p2.y - p1.y) * (p4.x - p3.x);
	}

	double Area(const trimesh::vec2& p1, const trimesh::vec2& p2, const trimesh::vec2& p3)
	{
		return Cross(p1, p2, p1, p3);
	}

	double fArea(const trimesh::vec2& p1, const trimesh::vec2& p2, const trimesh::vec2& p3)
	{
		return fabs(Area(p1, p2, p3));
	}

    trimesh::vec2 Inter(const trimesh::vec2& p1, const trimesh::vec2& p2, const trimesh::vec2& p3, const trimesh::vec2& p4)
	{
		double k = fArea(p1, p2, p3) / fArea(p1, p2, p4);
		return trimesh::vec2((p3.x + k * p4.x) / (1 + k), (p3.y + k * p4.y) / (1 + k));
	}

    void getHitData(openvdb::FloatGrid::Ptr gridptr, std::vector<trimesh::vec3>* supportPoints, std::vector<openvdb::math::Ray<double>>* rays,trimesh::TriMesh* omesh)
    {

		openvdb::OPENVDB_VERSION_NAME::tools::VolumeRayIntersector<openvdb::FloatGrid> aRayIntersector(*gridptr);
		for (openvdb::math::Ray<double>& aRay : *rays)
		{
			aRayIntersector.setIndexRay(aRay);
			std::list<openvdb::OPENVDB_VERSION_NAME::math::Ray<double>::TimeSpan> alist;
			aRayIntersector.hits(alist);

			openvdb::math::Vec3d cPointMin(omesh->bbox.max.x + 10, omesh->bbox.max.y + 10, omesh->bbox.max.z + 10);
			openvdb::math::Vec3d cPointMax(omesh->bbox.min.x - 10, omesh->bbox.min.y - 10, omesh->bbox.min.z - 10);
			for (openvdb::OPENVDB_VERSION_NAME::math::Ray<double>::TimeSpan atime : alist)
			{
				openvdb::math::Vec3d cPoint0 = aRayIntersector.getIndexPos(atime.t0);
				openvdb::math::Vec3d cPoint1 = aRayIntersector.getIndexPos(atime.t1);
				if (cPoint0[2] < cPointMin[2])
				{
					cPointMin = cPoint0;
				}
				if (cPoint1[2] < cPointMin[2])
				{
					cPointMin = cPoint1;
				}
				if (cPoint0[2] > cPointMax[2])
				{
					cPointMax = cPoint0;
				}
				if (cPoint1[2] > cPointMax[2])
				{
					cPointMax = cPoint1;
				}

			}
			supportPoints->push_back(trimesh::vec3(cPointMin[0], cPointMin[1], cPointMin[2]));
			supportPoints->push_back(trimesh::vec3(cPointMax[0], cPointMax[1], cPointMax[2]));
		}
    };

    std::vector<openvdb::math::Ray<double>>* generatorRay(trimesh::TriMesh* amesh, INNER_FILL_CONFIG fillConfig)
    {
        if (!fillConfig.enable)
        {
            return nullptr;
        }

        std::vector<openvdb::math::Ray<double>>* result = new std::vector<openvdb::math::Ray<double>>;
        std::vector<trimesh::vec2> Xlines;
        std::vector<trimesh::vec2> Ylines;

        amesh->need_bbox();
        float minZ = amesh->bbox.min.z;
        float minX = amesh->bbox.min.x;
        float maxX = amesh->bbox.max.x;
		float minY = amesh->bbox.min.y;
		float maxY = amesh->bbox.max.y;
        float gap = 20;

        int startX = minX / gap + 1;
        int endX = maxX / gap;
		int startY = minY / gap + 1;
		int endY = maxY / gap;

        if (startX*gap - minX<gap*0.2)
        {
            startX++;
        }
        if (maxX- endX * gap < gap * 0.2)
        {
            endX--;
        }

		if (startY * gap - minY < gap * 0.2)
		{
			startY++;
		}
		if (maxY - endY * gap < gap * 0.2)
		{
			endY--;
		}


        for (int n=startX;n<=endX;n++)
        {
            trimesh::vec2 pointStart(n *gap,minY);
            trimesh::vec2 pointEnd(n * gap, maxY);
            Xlines.push_back(pointStart);
            Xlines.push_back(pointEnd);
        }
		for (int n = startY; n <= endY; n++)
		{
			trimesh::vec2 pointStart(minX, n * gap);
			trimesh::vec2 pointEnd(maxX, n * gap);
            Ylines.push_back(pointStart);
            Ylines.push_back(pointEnd);
		}

        std::vector<trimesh::vec3> vctIntersection;
        for (int n=1;n<Xlines.size();n+=2)
        {
            //n-1,n
            for (int m=1;m<Ylines.size();m+=2)
            {
                trimesh::vec2 inters = Inter(Xlines[n - 1], Xlines[n], Ylines[m - 1], Ylines[m]);
                trimesh::point pointTemp(inters.x, inters.y,minZ);
                vctIntersection.push_back(pointTemp);
                
            }
        }

		using Vec3Type = openvdb::math::Vec3d;
		using Vec3T = Vec3Type;
        for (trimesh::vec3& apoint :vctIntersection)
        {
            Vec3Type eye = Vec3Type(apoint.at(0), apoint.at(1), apoint.at(2));
            Vec3Type direction = Vec3Type(0.0, 0.0, 1.0);
            openvdb::math::Ray<double> aRay(eye, direction);
            result->push_back(aRay);
        }

        return result;
    };


    static trimesh::TriMesh* _generate_interior(trimesh::TriMesh* mesh,
        double               min_thickness,
        double               voxel_scale,
        double               closing_dist,
        ccglobal::Tracer* tracer, double voxel_size, INNER_FILL_CONFIG fillConfig, std::vector<trimesh::vec3>* supportPoints)
    {
        //_scale(voxel_scale, imesh);

        double offset = voxel_scale * min_thickness;
        double D = voxel_scale * closing_dist;
        float  out_range = 0.03f * float(offset);
        float  in_range = 1.9f * float(offset + D);
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.0f);

        openvdb::FloatGrid::Ptr gridptr  = mesh_to_grid(mesh, {}, out_range, in_range, voxel_size);
        
        if (!gridptr) {
            if(tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.3f);

        if (closing_dist > .0) {
            gridptr = redistance_grid(*gridptr, -(offset + D), double(in_range));
        }
        else {
            D = -offset;
        }
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.7f);

        double iso_surface = D;
        double adaptivity = 0.;
        auto omesh = grid_to_mesh(*gridptr, iso_surface, adaptivity, false);

		//Моід
        omesh->need_bbox();
		std::vector<openvdb::math::Ray<double>>* rays = new std::vector<openvdb::math::Ray<double>>;
		rays = generatorRay(omesh, fillConfig);
        if (rays!=nullptr)
        {
            getHitData(gridptr, supportPoints,rays,omesh);
        }
        
        //_scale(1. / voxel_scale, omesh);
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(1.0f);
        
        return omesh;
    }

    trimesh::TriMesh* generateInterior(trimesh::TriMesh* mesh, std::vector<trimesh::vec3>* supportPoints,
        const HollowingParameter& parameter,ccglobal::Tracer* tracer)
    {
        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;

        // I can't figure out how to increase the grid resolution through openvdb
        // API so the model will be scaled up before conversion and the result
        // scaled down. Voxels have a unit size. If I set voxelSize smaller, it
        // scales the whole geometry down, and doesn't increase the number of
        // voxels.
        //
        // max 8x upscale, min is native voxel size
        // double voxel_scale = MIN_OVERSAMPL + (MAX_OVERSAMPL - MIN_OVERSAMPL) * parameter.quality;
        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* meshptr = _generate_interior(mesh, parameter.min_thickness, voxel_scale,
                parameter.closing_distance, tracer, parameter.voxel_size, parameter.fill_config, supportPoints);


        if (meshptr) {
            mmesh::reverseTriMesh(meshptr);
        }

        return meshptr;
    }

    static trimesh::TriMesh* _generate_boolcom(ovdbutil::TwoTrimesh* mesh,
        const int type, ccglobal::Tracer* tracer)
    {
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.0f);

        float offset = 1;
        float D = 0;
        float  out_range = 1.0f * float(offset);
        float  in_range = 1.f * float(offset + D);
        float voxel_size = 0.05;
        openvdb::FloatGrid::Ptr gridptr1 = mesh_to_grid(mesh->m1, {}, out_range, in_range, voxel_size);
        openvdb::FloatGrid::Ptr gridptr2 = mesh_to_grid(mesh->m2, {}, out_range, in_range, voxel_size);

        if (type==0)  openvdb::tools::csgIntersection(*gridptr1, *gridptr2);
        if (type == 1)  openvdb::tools::csgUnion(*gridptr1, *gridptr2);
        if (type == 2)  openvdb::tools::csgDifference(*gridptr1, *gridptr2);

        if (tracer)
            tracer->progress(0.4f);
        if (!gridptr1) {
            if (tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }
        if (!gridptr2) {
            if (tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }

        if (tracer)
            tracer->progress(0.7f);
        


        if (tracer)
            tracer->progress(1.0f);
        double iso_surface = 0.;
        double adaptivity = 0.;
        auto omesh = grid_to_mesh(*gridptr1, iso_surface, adaptivity, false);

        return omesh;
    }
    trimesh::TriMesh* generateBoolcom(ovdbutil::TwoTrimesh* mesh,
        const int type, ccglobal::Tracer* tracer)
    {
        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;

//        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* meshptr = _generate_boolcom(mesh, type,tracer);

        if (meshptr) {
            mmesh::reverseTriMesh(meshptr);
        }

        return meshptr;
    }

    void hollowMesh(trimesh::TriMesh* mesh,
        const HollowingParameter& parameter, ccglobal::Tracer* tracer)
    {

    }
}