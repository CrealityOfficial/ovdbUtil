#include "ovdbutil/booleanCom.h"

#include "util.h"
#include "mmesh/trimesh/trimeshutil.h"

#include "ccglobal/tracer.h"

namespace ovdbutil
{
    static trimesh::TriMesh* _generate_boolcom(ovdbutil::TwoTrimesh* mesh,
        const int type, const ovdbutil::BooleanParameter& param, ccglobal::Tracer* tracer)
    {
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.0f);

        float offset = 1;
        float D = 0;
        float  out_range = param.externalWidth;
        float  in_range = param.internalWidth;
        float voxel_size = 0.05;
        

        openvdb::initialize();

        // setup linear transform   
        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(0.25);

        // CUBE 
        std::vector<openvdb::math::Vec3s> cube_points;
        std::vector<openvdb::math::Coord::Vec3I> cube_faces;
        for (int i = 0; i < mesh->m1->vertices.size(); i++)
        {
            cube_points.push_back(openvdb::math::Vec3s(mesh->m1->vertices.at(i).x * 4, mesh->m1->vertices.at(i).y * 4, mesh->m1->vertices.at(i).z * 4));
        }
        for (int i = 0; i < mesh->m1->faces.size(); i++)
        {
            cube_faces.push_back(openvdb::math::Coord::Vec3I(mesh->m1->faces.at(i).x, mesh->m1->faces.at(i).y, mesh->m1->faces.at(i).z));
        }

        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_a(cube_points, cube_faces);
        //openvdb::FloatGrid::Ptr gridptr;
        openvdb::FloatGrid::Ptr gridptr1 = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_a, *xform);

        // spherer 
        std::vector<openvdb::math::Vec3s> s_points;
        std::vector<openvdb::math::Coord::Vec3I> s_faces;
        for (int i = 0; i < mesh->m2->vertices.size(); i++)
        {
            s_points.push_back(openvdb::math::Vec3s(mesh->m2->vertices.at(i).x * 4, mesh->m2->vertices.at(i).y * 4, mesh->m2->vertices.at(i).z * 4));
        }
        for (int i = 0; i < mesh->m2->faces.size(); i++)
        {
            s_faces.push_back(openvdb::math::Coord::Vec3I(mesh->m2->faces.at(i).x, mesh->m2->faces.at(i).y, mesh->m2->faces.at(i).z));
        }

        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_b(s_points, s_faces);
        //openvdb::FloatGrid::Ptr gridptr;
        openvdb::FloatGrid::Ptr gridptr2 = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_b, *xform);
        //openvdb::FloatGrid::Ptr gridptrout1 = mesh_to_grid(mesh->m1, {}, out_range, in_range, voxel_size);



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
        auto omesh = ovdbutil:: grid_to_mesh(*gridptr1, iso_surface, adaptivity, false);

        return omesh;
    }
    trimesh::TriMesh* generateBoolcom(ovdbutil::TwoTrimesh* mesh,
        const int type,const ovdbutil::BooleanParameter& param, ccglobal::Tracer* tracer)
    {
        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;

//        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* meshptr = _generate_boolcom(mesh, type,param,tracer);

        if (meshptr) {
            mmesh::reverseTriMesh(meshptr);
        }

        return meshptr;
    }


}