#ifndef OVDBUTIL_HOLLOWING_1650957593077_H
#define OVDBUTIL_HOLLOWING_1650957593077_H
#include "ovdbutil/interface.h"

namespace trimesh
{
    class TriMesh;
}

namespace ccglobal
{
    class Tracer;
}

namespace ovdbutil
{
    struct HollowingParameter
    {
        double min_thickness = 1.0;
        double quality = 0.5;
        double closing_distance = 0.0;
        double voxel_size_inout_range = 1.0;
        double voxel_size = 1.0;
    };

    struct TwoTrimesh
    {
        trimesh::TriMesh *m1;
        trimesh::TriMesh *m2;
    };

    OVDBUTIL_API trimesh::TriMesh* generateBoolcom(ovdbutil::TwoTrimesh* mesh, const int type, ccglobal::Tracer* tracer=nullptr);
    OVDBUTIL_API trimesh::TriMesh* generateInterior(trimesh::TriMesh* mesh,
        const HollowingParameter& = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    OVDBUTIL_API void hollowMesh(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);
}

#endif // OVDBUTIL_HOLLOWING_1650957593077_H