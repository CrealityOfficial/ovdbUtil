#ifndef OVDBUTIL_BOOLEAN_1650957593077_H
#define OVDBUTIL_BOOLEAN_1650957593077_H
#include "ovdbutil/interface.h"
#include "mmesh/trimesh/trimeshutil.h"
#include <vector>
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
    OVDBUTIL_API trimesh::TriMesh* generateDilationcom(trimesh::TriMesh* mesh, const int param, ccglobal::Tracer* tracer = nullptr);
    OVDBUTIL_API trimesh::TriMesh* generateErosioncom(trimesh::TriMesh* mesh, const int param, ccglobal::Tracer* tracer = nullptr);
    OVDBUTIL_API trimesh::TriMesh* voroniaBoolcom(std::vector<trimesh::TriMesh>* mesh, const int type, const ovdbutil::BooleanParameter& param, ccglobal::Tracer* tracer=nullptr);
    OVDBUTIL_API trimesh::TriMesh* voroniaBoolcomInputLines(std::vector<edge> &edges);
}

#endif // OVDBUTIL_HOLLOWING_1650957593077_H