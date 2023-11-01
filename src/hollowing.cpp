#include "ovdbutil/hollowing.h"

#include "util.h"

#include "mesh/trimeshutil.h"
#include "mesh/createcylinder.h"
#include "mesh/ray.h"
#include "mesh/adjacentoctree.h"

#include "topomesh/interface/subdivision.h"
#include "trimesh2/TriMesh_algo.h"
#include <openvdb/math/Vec3.h>
#include <openvdb/math/Coord.h>
#include <openvdb/tools/LevelSetSphere.h>

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
       // float  out_range = 0.01f * float(offset);
       float  in_range = 1.9f * float(offset + D);
        //float  in_range = 2.8f * float(offset + D);
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
        //openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(/*background value=*/2.0);

        //bug for param flags: hollow Optimization 
       openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(0.333);
        //std::cout << "trans size :" << transform->voxelSize() << "\n";
        //transform.postScale(0.5);       
       // std::cout << "after trans size :" << transform->voxelSize() << "\n";

        std::vector<openvdb::math::Vec3s> cube_points;
        std::vector<openvdb::math::Coord::Vec3I> cube_faces;
        for (int i = 0; i < mesh->vertices.size(); i++)
        {
            cube_points.push_back(openvdb::math::Vec3s(mesh->vertices.at(i).x*3 , mesh->vertices.at(i).y*3, mesh->vertices.at(i).z*3));
        }
        for (int i = 0; i < mesh->faces.size(); i++)
        {
            cube_faces.push_back(openvdb::math::Coord::Vec3I(mesh->faces.at(i).x, mesh->faces.at(i).y, mesh->faces.at(i).z));
        }
        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_a(cube_points, cube_faces);
        openvdb::FloatGrid::Ptr gridptr1 = mesh_to_grid(mesh, {}, out_range, in_range, voxel_size, 0xE, tracer);
        //openvdb::FloatGrid::Ptr gridptr1 = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_a, *transform, out_range, in_range, voxel_size, 0xE, tracer);
        
        //openvdb::FloatGrid::Ptr gridptr =
        //    openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
        //        /*radius=*/10.0, /*center=*/openvdb::Vec3f(0, 0, 0),
        //        /*voxel size=*/0.5, /*width=*/3.0);   
        //openvdb::Vec3f c = openvdb::Vec3f(0, 0, 0);
        //openvdb::FloatGrid::Ptr grid =
        //    openvdb::FloatGrid::create(/*background value=*/2.0);
        //using ValueT = typename openvdb::FloatGrid::ValueType;
        //const ValueT outside = grid->background();
        //const ValueT inside = -outside;
        //int padding = int(openvdb::math::RoundUp(openvdb::math::Abs(outside)));
        //std::cout << "background : " << grid->background() << "\n";
        //std::cout << "padding : " << padding << "\n";
        //int dim = int(10.f + padding);
        //typename openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
        //openvdb::Coord ijk;
        //int& i = ijk[0], & j = ijk[1], & k = ijk[2];
        //for (i = c[0] - dim; i < c[0] + dim; ++i) {
        //    const float x2 = openvdb::math::Pow2(i - c[0]);
        //    for (j = c[1] - dim; j < c[1] + dim; ++j) {
        //        const float x2y2 = openvdb::math::Pow2(j - c[1]) + x2;
        //        for (k = c[2] - dim; k < c[2] + dim; ++k) {
        //            // The distance from the sphere surface in voxels
        //            const float dist = openvdb::math::Sqrt(x2y2
        //                + openvdb::math::Pow2(k - c[2])) - 10.f;
        //            // Convert the floating-point distance to the grid's value type.
        //            ValueT val = ValueT(dist);
        //            // Only insert distances that are smaller in magnitude than
        //            // the background value.
        //            if (val < inside || outside < val) continue;
        //            // Set the distance for voxel (i,j,k).
        //            accessor.setValue(ijk, val);
        //        }
        //    }
        //}
        //openvdb::tools::signedFloodFill(grid->tree());
       /* for (openvdb::FloatGrid::ValueOnIter iter = gridptr->beginValueOn(); iter; ++iter) {
            float dist = iter.getValue();
            iter.setValue((outside - dist) / width);
        }*/
        //gridptr->insertMeta("radius", openvdb::FloatMetadata(50.0));
        //std::cout << "Grid size :" << gridptr->voxelSize()<<" background : " <<gridptr->background();
       // gridptr->transformPtr()->postScale(0.5f);
       /* for (openvdb::FloatGrid::ValueOnIter iter = gridptr->beginValueOn(); iter; ++iter) {
            float dist = iter.getValue();
            iter.setValue(dist/2.0f);
        }*/
       // std::cout << "Grid size :" << gridptr->voxelSize() << " background : " << gridptr->background();
       /* for (openvdb::FloatGrid::ValueOnCIter iter = gridptr->cbeginValueOn(); iter; ++iter) {
            std::cout << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
        }*/
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


        if (!gridptr1) {
            if(tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (closing_dist > .0) {
            gridptr1 = redistance_grid(*gridptr1, -(offset + D), double(in_range));
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
        //double iso_surface = 0;
        double adaptivity = 0.;
        auto omesh = grid_to_mesh(*gridptr1, iso_surface, adaptivity, false);
        
        //_scale(1. / voxel_scale, omesh);
        
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.95f);
        
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
        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* meshptr = _generate_interior(mesh, parameter.min_thickness, voxel_scale,
                parameter.closing_distance, tracer, parameter.voxel_size);


        if (meshptr) {
            mmesh::reverseTriMesh(meshptr);
        }

        return meshptr;
    }

    trimesh::TriMesh* hollowPrecisionMeshAndFill(trimesh::TriMesh* mesh,
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
        trimesh::TriMesh* returnmesh = hollowMesh(mesh, parameter, tracer);
        
        trimesh::TriMesh* outMesh = new trimesh::TriMesh();
        std::vector<trimesh::TriMesh*> meshtotalV;
        meshtotalV.push_back(mesh);
        if (returnmesh != nullptr && returnmesh->vertices.size() > 0)
        {
            meshtotalV.push_back(returnmesh);
            if (1)//
            {
                trimesh::vec3 normal(0.5, 0.0, 1.0);
                std::vector<trimesh::TriMesh*> vctMesh = generateInfill(returnmesh, normal, parameter);
                for (trimesh::TriMesh* cyMesh : vctMesh)
                {
                    meshtotalV.push_back(cyMesh);
                }
            }
            if (2)//
            {
                trimesh::vec3 normal(-0.5, 0.0, 1.0);
                std::vector<trimesh::TriMesh*> vctMesh = generateInfill(returnmesh, normal, parameter);
                for (trimesh::TriMesh* cyMesh : vctMesh)
                {
                    meshtotalV.push_back(cyMesh);
                }
            }
            if (3)//
            {
                trimesh::vec3 normal(0.0, 0.5, 1.0);
                std::vector<trimesh::TriMesh*> vctMesh = generateInfill(returnmesh, normal, parameter);
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
       // mesh->write("origin.ply");
       // trimesh::TriMesh* newmesh = new trimesh::TriMesh();
       // mesh->need_normals();
       // mesh->need_neighbors();
       // mesh->need_adjacentfaces();
       // for (int vi = 0; vi < mesh->vertices.size(); vi++)
       // {           
       //     if (vi == 67511)
       //         std::cout << "";
       //     trimesh::point n(0, 0, 0);
       //     for (int fi = 0; fi < mesh->adjacentfaces[vi].size(); fi++)
       //     {
       //         int f = mesh->adjacentfaces[vi][fi];
       //         trimesh::point nt = (mesh->vertices[mesh->faces[f][1]] - mesh->vertices[mesh->faces[f][0]]) % (mesh->vertices[mesh->faces[f][2]] - mesh->vertices[mesh->faces[f][0]]);
       //         trimesh::normalize(nt);
       //         n += nt;
       //     }
       //     n /= mesh->adjacentfaces[vi].size();
       //     trimesh::point t1(mesh->vertices[vi]-n);
       //     /*trimesh::point c(0, 0, 0);
       //     for (int vv = 0; vv < mesh->neighbors[vi].size(); vv++)
       //     {
       //         c += mesh->normals[mesh->neighbors[vi][vv]];
       //     }
       //     c /= mesh->neighbors[vi].size();
       //     trimesh::point t2(mesh->vertices[vi] - c);*/
       //     newmesh->vertices.push_back(t1);
       // }
       // for (int fi = 0; fi < mesh->faces.size(); fi++)
       // {
       //     newmesh->faces.push_back(mesh->faces[fi]);
       // }
       // trimesh::TriMesh* pointmesh = new trimesh::TriMesh();
       // newmesh->need_normals();
       // newmesh->need_adjacentfaces();
       // std::vector<bool> del(newmesh->vertices.size(), false);
       // for (int vi = 0; vi < newmesh->vertices.size(); vi++)
       // {
       //     trimesh::point n1 = newmesh->normals[vi];
       //     trimesh::point n2 = mesh->normals[vi];
       //     float arc = n1 ^ n2;
       //     arc = arc >= 1.f ? 1.f : arc;
       //     arc = arc <= -1.f ? -1.f : arc;
       //     float ang = std::acos(arc) * 180 / M_PI;
       //     if (ang > 45||ang<-135)
       //     {
       //         pointmesh->vertices.push_back(newmesh->vertices[vi]);
       //         del[vi] = true;
       //       /*  for (int fi = 0; fi < newmesh->adjacentfaces[vi].size(); fi++)
       //         {
       //             delf[newmesh->adjacentfaces[vi][fi]] = true;
       //         }*/
       //     }
       // }
       // trimesh::TriMesh* savemesh = new trimesh::TriMesh();
       // for (int fi = 0; fi < newmesh->faces.size(); fi++)
       // {
       //     int v0 = newmesh->faces[fi][0];
       //     int v1 = newmesh->faces[fi][1];
       //     int v2 = newmesh->faces[fi][2];
       //     if (del[v0] || del[v1] || del[v2])
       //         continue;
       //     savemesh->vertices.push_back(newmesh->vertices[v0]);
       //     savemesh->vertices.push_back(newmesh->vertices[v1]);
       //     savemesh->vertices.push_back(newmesh->vertices[v2]);
       //     savemesh->faces.push_back(trimesh::ivec3(savemesh->vertices.size()-3, savemesh->vertices.size() - 2, savemesh->vertices.size() - 1));
       // }
       // //trimesh::remove_faces(newmesh,delf);
       //// trimesh::remove_vertices(newmesh,del);
       //// trimesh::remove_unused_vertices(newmesh);
       ///* for (int fi = 0; fi < mesh->faces.size(); fi++)
       // {
       //     trimesh::point n = (mesh->vertices[mesh->faces[fi][1]] - mesh->vertices[mesh->faces[fi][0]])
       //         % (mesh->vertices[mesh->faces[fi][2]] - mesh->vertices[mesh->faces[fi][0]]);
       //     trimesh::normalize(n);
       //     trimesh::point v0 = mesh->vertices[mesh->faces[fi][0]] - n;
       //     trimesh::point v1 = mesh->vertices[mesh->faces[fi][1]] - n;
       //     trimesh::point v2 = mesh->vertices[mesh->faces[fi][2]] - n;
       //     newmesh->vertices.push_back(v0);
       //     newmesh->vertices.push_back(v1);
       //     newmesh->vertices.push_back(v2);
       //     newmesh->faces.push_back(trimesh::ivec3(fi*3,fi*3+1,fi*3+2));
       // }*/
       // savemesh->write("savemesh.ply");
       // pointmesh->write("pointmesh.ply");
       // newmesh->write("newmesh.ply");
        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;
        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* hollowMesh = _generate_interior(mesh, parameter.min_thickness, voxel_scale,
            parameter.closing_distance, tracer, parameter.voxel_size);
        hollowMesh->write("hollowstep0.ply");
      /*  std::vector<int> faceindex(hollowMesh->faces.size());
        for (int i = 0; i < hollowMesh->faces.size(); i++)
            faceindex[i] = i;
        std::vector<int> resultfaces;
        for (int i = 0; i < 2; i++)
        {
            topomesh::loopSubdivision(hollowMesh, faceindex, resultfaces);
            faceindex.swap(resultfaces);
            resultfaces.clear();
        }
        hollowMesh->need_normals();
        std::random_device rd;
        std::mt19937 engine(rd());
        std::uniform_real_distribution<double> spans(0.0, 0.2);
        for (int vi = 0; vi < hollowMesh->vertices.size(); vi++)
        {
            hollowMesh->vertices[vi] -= spans(engine)*trimesh::normalized(hollowMesh->normals[vi]);
        }
        hollowMesh->need_neighbors();
        for (int vi = 0; vi < hollowMesh->vertices.size(); vi++)
        {
            trimesh::point all(0, 0, 0);
            for (int vv = 0; vv < hollowMesh->neighbors[vi].size(); vv++)
                all += hollowMesh->vertices[hollowMesh->neighbors[vi][vv]];
            all /= hollowMesh->neighbors[vi].size();
            hollowMesh->vertices[vi] = all;
        }*/

       
       /* for (int vi = 0; vi < hollowMesh->vertices.size(); vi++)
        {
            trimesh::point c(0, 0, 0);
            for (int vv = 0; vv < hollowMesh->neighbors[vi].size(); vv++)
            {
                c += hollowMesh->vertices[hollowMesh->neighbors[vi][vv]];
            }
            c /= hollowMesh->neighbors[vi].size();
            hollowMesh->vertices[vi] = hollowMesh->vertices[vi] - c;
        }*/
       
       
        //抽壳后填充0
   		trimesh::TriMesh* outMesh = new trimesh::TriMesh();
		std::vector<trimesh::TriMesh*> meshtotalV;
		meshtotalV.push_back(mesh);
        if (hollowMesh != nullptr && hollowMesh->vertices.size() > 0)
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
        hollowMesh->write("hollowstep1.ply");
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

    trimesh::TriMesh* hollowMesh(trimesh::TriMesh* mesh,
        const HollowingParameter& parameter, ccglobal::Tracer* tracer)
    {
        float min_thickness = parameter.min_thickness * parameter.precision*1.0f;
        double voxel_scale = parameter.voxel_size_inout_range;
        double offset = voxel_scale * min_thickness;
        double D = voxel_scale * parameter.closing_distance;
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
        openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(1.0f/ (parameter.precision * 1.0f));
        std::vector<openvdb::math::Vec3s> cube_points;
        std::vector<openvdb::math::Coord::Vec3I> cube_faces;
        for (int i = 0; i < mesh->vertices.size(); i++)
        {
            cube_points.push_back(openvdb::math::Vec3s(mesh->vertices.at(i).x * parameter.precision, mesh->vertices.at(i).y * parameter.precision,
                mesh->vertices.at(i).z * parameter.precision));
        }
        for (int i = 0; i < mesh->faces.size(); i++)
        {
            cube_faces.push_back(openvdb::math::Coord::Vec3I(mesh->faces.at(i).x, mesh->faces.at(i).y, mesh->faces.at(i).z));
        }
        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_a(cube_points, cube_faces);
        openvdb::FloatGrid::Ptr gridptr = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_a, *transform, out_range, in_range, parameter.voxel_size, 0xE, tracer);
        if (!gridptr) {
            if (tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }

        if (tracer && tracer->interrupt())
            return nullptr;

        if (parameter.closing_distance > .0) {
            gridptr = redistance_grid(*gridptr, -(offset + D), double(in_range));
        }
        else {
            D = -offset;
        }

        double iso_surface = D;
        double adaptivity = 0.;
        trimesh::TriMesh* hollowMesh = grid_to_mesh(*gridptr, iso_surface, adaptivity, false);

        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.95f);
      
        return hollowMesh;
    }
}
