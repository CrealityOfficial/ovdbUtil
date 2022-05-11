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
        openvdb::FloatGrid::Ptr gridptr1 = ovdbutil::mesh_to_grid(mesh->m1, {}, out_range, in_range, voxel_size);
        //openvdb::FloatGrid::Ptr gridptrout1 = mesh_to_grid(mesh->m1, {}, out_range, in_range, voxel_size);
        openvdb::FloatGrid::Ptr gridptr2 = ovdbutil::mesh_to_grid(mesh->m2, {}, out_range, in_range, voxel_size);
        //openvdb::FloatGrid::Ptr gridptrout2 = mesh_to_grid(mesh->m2, {}, out_range, in_range, voxel_size);
        //// Get the source and target grids' index space to world space transforms.
        // const openvdb::math::Transform
        //     & sourceXform = gridptr1->transform(),
        //     & targetXform = gridptrout1->transform();
        // // Compute a source grid to target grid transform.
        // // (For this example, we assume that both grids' transforms are linear,
        // // so that they can be represented as 4 x 4 matrices.)
        // openvdb::Mat4R xform =
        //     sourceXform.baseMap()->getAffineMap()->getMat4() *
        //     targetXform.baseMap()->getAffineMap()->getMat4().inverse();
        // //for (int i=0; i< 16; i++)
        // //    xform.asPointer()[i]= xform.asPointer()[i]*0.01;
        // // Create the transformer.
        // openvdb::tools::GridTransformer transformer(xform);
        // // Resample using nearest-neighbor interpolation.
        // transformer.transformGrid<openvdb::tools::QuadraticSampler, openvdb::FloatGrid>(
        //     *gridptr1, *gridptrout1);
        // // Prune the target tree for optimal sparsity.
        // gridptrout1->tree().prune();


        // // Get the source and target grids' index space to world space transforms.
        // const openvdb::math::Transform
        //     & sourceXform2 = gridptr2->transform(),
        //     & targetXform2 = gridptrout2->transform();
        // // Compute a source grid to target grid transform.
        // // (For this example, we assume that both grids' transforms are linear,
        // // so that they can be represented as 4 x 4 matrices.)
        // openvdb::Mat4R xform2 =
        //     sourceXform2.baseMap()->getAffineMap()->getMat4() *
        //     targetXform2.baseMap()->getAffineMap()->getMat4().inverse();
        // //for (int i=0; i< 16; i++)
        // //    xform.asPointer()[i]= xform.asPointer()[i]*0.01;
        // // Create the transformer.
        // openvdb::tools::GridTransformer transformer2(xform2);
        // // Resample using nearest-neighbor interpolation.
        // transformer2.transformGrid<openvdb::tools::QuadraticSampler, openvdb::FloatGrid>(
        //     *gridptr2, *gridptrout2);
        // // Prune the target tree for optimal sparsity.
        // gridptrout2->tree().prune();


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