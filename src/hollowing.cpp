#include "ovdbutil/hollowing.h"
#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetRebuild.h>

#include "trimesh2/TriMesh.h"
#include "mmesh/trimesh/trimeshutil.h"

#include "ccglobal/tracer.h"

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
        int                 flags = 0)
    {
        openvdb::initialize();

        openvdb::FloatGrid::Ptr grid = openvdb::tools::meshToVolume<openvdb::FloatGrid>(
                TriangleMeshDataAdapter{ *mesh }, tr, exteriorBandWidth,
                interiorBandWidth, flags);

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

    static trimesh::TriMesh* _generate_interior(trimesh::TriMesh* mesh,
        double               min_thickness,
        double               voxel_scale,
        double               closing_dist,
        ccglobal::Tracer* tracer)
    {
        //_scale(voxel_scale, imesh);

        double offset = voxel_scale * min_thickness;
        double D = voxel_scale * closing_dist;
        float  out_range = 0.1f * float(offset);
        float  in_range = 1.1f * float(offset + D);
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.0f);

        openvdb::FloatGrid::Ptr gridptr  = mesh_to_grid(mesh, {}, out_range, in_range);
        
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
        
        //_scale(1. / voxel_scale, omesh);
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(1.0f);
        
        return omesh;
    }

    trimesh::TriMesh* generateInterior(trimesh::TriMesh* mesh,
        const HollowingParameter& parameter, ccglobal::Tracer* tracer)
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
        double voxel_scale = 1.0;
        trimesh::TriMesh* meshptr = _generate_interior(mesh, parameter.min_thickness, voxel_scale,
                parameter.closing_distance, tracer);

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