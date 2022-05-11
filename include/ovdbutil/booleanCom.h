#ifndef OVDBUTIL_BOOLEAN_1650957593077_H
#define OVDBUTIL_BOOLEAN_1650957593077_H
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
    struct TwoTrimesh
    {
        trimesh::TriMesh *m1;
        trimesh::TriMesh *m2;
    };

    struct BooleanParameter
    {
        double externalWidth;
        double internalWidth;
    };

    OVDBUTIL_API trimesh::TriMesh* generateBoolcom(ovdbutil::TwoTrimesh* mesh, const int type, const BooleanParameter& param, ccglobal::Tracer* tracer=nullptr);
}

#endif // OVDBUTIL_HOLLOWING_1650957593077_H