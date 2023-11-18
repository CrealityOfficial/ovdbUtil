#ifndef OVDBUTIL_GRID_HELPER_UTIL_H
#define OVDBUTIL_GRID_HELPER_UTIL_H
#include "dataadapter.h"
#include <openvdb/openvdb.h>

namespace ovdbutil
{
    typedef openvdb::FloatGrid::Ptr FloatGridPtr;

    class TracerInterrupter;
    FloatGridPtr mesh_to_grid(
        TracerInterrupter& tracer,
        trimesh::TriMesh* mesh,
        const openvdb::math::Transform& tr,
        float               exteriorBandWidth,
        float               interiorBandWidth,
        double              voxel_size,
        int                 flags = 0
    );

    trimesh::TriMesh* grid_to_mesh(FloatGridPtr grid,
        double                    isovalue,
        double                    adaptivity,
        bool                      relaxDisorientedTriangles);

    Contour3D _volumeToMesh(FloatGridPtr grid,
        double      isovalue,
        double      adaptivity,
        bool        relaxDisorientedTriangles);
}

#endif // OVDBUTIL_GRID_HELPER_UTIL_H