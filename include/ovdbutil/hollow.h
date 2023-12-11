#ifndef OVDBUTIL_HOLLOW_1650957593077_H
#define OVDBUTIL_HOLLOW_1650957593077_H
#include "ovdbutil/interface.h"
#include <vector>

namespace ovdbutil
{
    enum class OutPrecesion {
        COARSE = 0x1,
        NORMAL = 0x2,
        ELABORATE = 0x3,
        EXTREME = 0x4
    };

    enum class HollowStyle
    {
        hs_none,
        hs_infill_grid,
        hs_self_support
    };
   
    struct HollowParameter
    {
        double thickness = 2.0;
        OutPrecesion precision = OutPrecesion::NORMAL;

        bool remain_main_shell = true;
        bool filter_shell = false;
        double filter_tiny_shell = 5.0;  // mm3
        HollowStyle style = HollowStyle::hs_none;
        double fill_ratio = 20;
    };

    struct ShellParameter
    {
        double thickness = 2.0;
        OutPrecesion precision = OutPrecesion::NORMAL;

        std::vector<int> faces; // removed faces
    };

    OVDBUTIL_API trimesh::TriMesh* hollowMesh(trimesh::TriMesh* mesh,
        const HollowParameter& = HollowParameter(), ccglobal::Tracer* tracer = nullptr);

    OVDBUTIL_API trimesh::TriMesh* shellMesh(trimesh::TriMesh* mesh, const ShellParameter& = ShellParameter(), ccglobal::Tracer* tracer = nullptr);
}

#endif // OVDBUTIL_HOLLOW_1650957593077_H