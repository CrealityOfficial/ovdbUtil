#include "ovdbutil/hollowing.h"

#include "util.h"
#include "mmesh/trimesh/trimeshutil.h"
#include "ccglobal/tracer.h"

namespace ovdbutil
{
    openvdb::FloatGrid::Ptr redistance_grid(const openvdb::FloatGrid& grid, double iso, double er = 3.0, double ir = 3.0)
    {
        return openvdb::tools::levelSetRebuild(grid, float(iso), float(er), float(ir));
    }

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
            tracer->progress(0.0f);

        openvdb::FloatGrid::Ptr gridptr  = mesh_to_grid(mesh, {}, out_range, in_range, voxel_size);
        


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
        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* meshptr = _generate_interior(mesh, parameter.min_thickness, voxel_scale,
                parameter.closing_distance, tracer, parameter.voxel_size);

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