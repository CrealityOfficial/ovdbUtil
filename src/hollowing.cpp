#include "ovdbutil/hollowing.h"

#include "util.h"
#include "mmesh/trimesh/trimeshutil.h"
#include "ccglobal/tracer.h"
#include <openvdb/tools/RayIntersector.h>
#include "mmesh/create/createcylinder.h"
#include "mmesh/camera/ray.h"
#include "mmesh/util/adjacentoctree.h"
#include <openvdb/math/Vec3.h>
#include <openvdb/math/Coord.h>
namespace ovdbutil
{
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

    std::vector<mmesh::Ray>* generatorRay(trimesh::TriMesh* amesh, INNER_FILL_CONFIG fillConfig, const trimesh::vec3& normal)
    {

        std::vector<mmesh::Ray>* result = new std::vector<mmesh::Ray>;
        std::vector<trimesh::vec2> Xlines;
        std::vector<trimesh::vec2> Ylines;

        amesh->need_bbox();
        float gap = fillConfig.fillRadius / fillConfig.fillratio;
        float minZ = amesh->bbox.min.z;
        float minX = amesh->bbox.min.x;
        float maxX = amesh->bbox.max.x;
		int startX = minX / gap + 1;
		int endX = maxX / gap;
		float minY = amesh->bbox.min.y;
		float maxY = amesh->bbox.max.y;
		int startY = minY / gap + 1;
		int endY = maxY / gap;

        //排除超出边界的交点
        if (startX* gap - minX< fillConfig.fillRadius*2)
        {
            startX++;
        }
        if (maxX- endX * gap < fillConfig.fillRadius * 2)
        {
            endX--;
        }
		if (startY * gap - minY < fillConfig.fillRadius * 2)
		{
			startY++;
		}
		if (maxY - endY * gap < fillConfig.fillRadius * 2)
		{
			endY--;
		}


        for (int n=startX;n<=endX;n++)
        {
            trimesh::vec2 pointStart(n * gap,minY);
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

        for (trimesh::vec3& apoint :vctIntersection)
        {
            mmesh::Ray aRay(apoint,normal);
            result->push_back(aRay);
        }
        return result;
    };


    static trimesh::TriMesh* _generate_interior(trimesh::TriMesh* mesh,
        double               min_thickness,
        double               voxel_scale,
        double               closing_dist,
        ccglobal::Tracer* tracer, double voxel_size)
    {
        //_scale(voxel_scale, imesh);

        double offset = voxel_scale * min_thickness;
        double D = voxel_scale * closing_dist;
        float  out_range = 0.03f * float(offset);
        float  in_range = 1.9f * float(offset + D);
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
        {
            tracer->progress(0.15f);
            if (tracer->interrupt())
            {
                return nullptr;
            }
        }


        //bug for param flags: hollow Optimization 
        openvdb::FloatGrid::Ptr gridptr  = mesh_to_grid(mesh, {}, out_range, in_range, voxel_size,0xE, tracer);
        
        //openvdb::FloatGrid::Ptr gridptrout1 = mesh_to_grid(mesh, {}, out_range, in_range, voxel_size);
        //// Get the source and target grids' index space to world space transforms.
        //const openvdb::math::Transform
        //    & sourceXform = gridptr->transform(),
        //    & targetXform = gridptrout1->transform();
        //// Compute a source grid to target grid transform.
        //// (For this example, we assume that both grids' transforms are linear,
        //// so that they can be represented as 4 x 4 matrices.)
        //openvdb::Mat4R xform =
        //    sourceXform.createLinearTransform(40.0).get()->baseMap()->getAffineMap()->getMat4() *
        //    targetXform.createLinearTransform(40.0).get()->baseMap()->getAffineMap()->getMat4().inverse();
        //for (int i=0; i< 16; i++)
        //    xform.asPointer()[i]= xform.asPointer()[i]*10;
        //// Create the transformer.
        //openvdb::tools::GridTransformer transformer(xform);
        //// Resample using nearest-neighbor interpolation.
        //transformer.transformGrid<openvdb::tools::QuadraticSampler, openvdb::FloatGrid>(
        //    *gridptr, *gridptrout1);
        //// Prune the target tree for optimal sparsity.
        //gridptrout1->tree().prune();


        if (!gridptr) {
            if(tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (closing_dist > .0) {
            gridptr = redistance_grid(*gridptr, -(offset + D), double(in_range));
        }
        else {
            D = -offset;
        }   

        if (tracer)
        {
            tracer->progress(0.8f);
            if (tracer->interrupt())
            {
                return nullptr;
            }
        }

        double iso_surface = D;
        double adaptivity = 0.;
        auto omesh = grid_to_mesh(*gridptr, iso_surface, adaptivity, false);
        
        //_scale(1. / voxel_scale, omesh);
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.95f);
        
        return omesh;
    }

    trimesh::TriMesh* generateInterior(trimesh::TriMesh* mesh, std::vector<trimesh::vec3>* supportPoints,
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
        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* meshptr = _generate_interior(mesh, parameter.min_thickness, voxel_scale,
                parameter.closing_distance, tracer, parameter.voxel_size);


        if (meshptr) {
            mmesh::reverseTriMesh(meshptr);
        }

        return meshptr;
    }

    OVDBUTIL_API trimesh::TriMesh* hollowMeshAndFill(trimesh::TriMesh* mesh,
        const HollowingParameter & parameter, ccglobal::Tracer* tracer)
    {
        if (tracer)
        {
            tracer->progress(0.1f);
            if (tracer->interrupt())
            {
                return nullptr;
            }
        }

        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;
        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* hollowMesh = _generate_interior(mesh, parameter.min_thickness, voxel_scale,
            parameter.closing_distance, tracer, parameter.voxel_size);

        //抽壳后填充0
		trimesh::TriMesh* outMesh = new trimesh::TriMesh();
		std::vector<trimesh::TriMesh*> meshtotalV;
		meshtotalV.push_back(mesh);
        if (hollowMesh) 
        {
            meshtotalV.push_back(hollowMesh);
            if (1)//
            {
				trimesh::vec3 normal(0.5, 0.0, 1.0);
				std::vector<trimesh::TriMesh*> vctMesh = generateInfill(hollowMesh, normal, parameter);
				for (trimesh::TriMesh* cyMesh : vctMesh)
				{
					meshtotalV.push_back(cyMesh);
				}
            }
			if (2)//
			{
				trimesh::vec3 normal(-0.5, 0.0, 1.0);
				std::vector<trimesh::TriMesh*> vctMesh = generateInfill(hollowMesh, normal, parameter);
				for (trimesh::TriMesh* cyMesh : vctMesh)
				{
					meshtotalV.push_back(cyMesh);
				}
			}
			if (3)//
			{
				trimesh::vec3 normal(0.0, 0.5, 1.0);
				std::vector<trimesh::TriMesh*> vctMesh = generateInfill(hollowMesh, normal, parameter);
				for (trimesh::TriMesh* cyMesh : vctMesh)
				{
					meshtotalV.push_back(cyMesh);
				}
			}
        }

        mmesh::mergeTriMesh(outMesh, meshtotalV);

        if (tracer)
        {
            tracer->progress(1.0f);
            if (tracer->interrupt())
            {
                return nullptr;
            }
        }
        return outMesh;
    }

    OVDBUTIL_API std::vector<trimesh::TriMesh*> generateInfill(trimesh::TriMesh* mesh, const trimesh::vec3& normal, const HollowingParameter& param, ccglobal::Tracer* tracer /*= nullptr*/)
	{
        std::vector<trimesh::TriMesh*> vctMesh;
        if (!param.fill_config.enable)
        {
            return vctMesh;
        }

        std::vector<mmesh::Ray>* rays = generatorRay(mesh, param.fill_config,normal);
        ModelAdjacentOctree aModelAdj(mesh);
        for (mmesh::Ray& aray :*rays)
        {
            std::vector<std::pair<int, trimesh::dvec3>> faceIdIntersect;
            aModelAdj.lineCollide((trimesh::dvec3)aray.start, (trimesh::dvec3)trimesh::normalized(aray.dir), faceIdIntersect);

            //过滤相同的点
            std::vector<trimesh::dvec3> vctIntersect;
            for (std::pair<int,trimesh::dvec3>& apair:faceIdIntersect)
            {
                static double preZ = apair.second.z - 1;
                if (preZ == apair.second.z)
                {
                    continue;
                }
                vctIntersect.push_back(apair.second);
                preZ = apair.second.z;
            }

            for (int n=1;n< vctIntersect.size();n+=2)
            {
                const trimesh::vec3& startPoint = (trimesh::vec3)vctIntersect.at(n-1);
                const trimesh::vec3& endPoint = (trimesh::vec3)vctIntersect.at(n);
				float height = trimesh::distance(startPoint, endPoint) + param.min_thickness * 0.5;//param.min_thickness*0.25填充柱子伸出一点，更好的附着
				trimesh::vec3 centerPoint((endPoint.x + startPoint.x) * 0.5, (endPoint.y + startPoint.y) * 0.5, (endPoint.z + startPoint.z) * 0.5);
				trimesh::TriMesh* cylinderMesh = mmesh::createSoupCylinder(10, param.fill_config.fillRadius, height, centerPoint, aray.dir);
				vctMesh.push_back(cylinderMesh);
            }
        }
        return vctMesh;
	}

	void hollowMesh(trimesh::TriMesh* mesh,
        const HollowingParameter& parameter, ccglobal::Tracer* tracer)
    {

    }
}
