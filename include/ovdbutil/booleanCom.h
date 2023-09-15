#ifndef OVDBUTIL_BOOLEAN_1650957593077_H
#define OVDBUTIL_BOOLEAN_1650957593077_H
#include "ovdbutil/interface.h"
#include "trimesh2/TriMesh.h"
#include "ccglobal/tracer.h"
#include <vector>

namespace ovdbutil
{
    struct edge {
        trimesh::vec3 p0;
        trimesh::vec3 p1;
    };

	/// @brief TwoTrimesh
    struct TwoTrimesh
    {
        trimesh::TriMesh *m1;
        trimesh::TriMesh *m2;
    };
	/// @brief BooleanParameter

    struct BooleanParameter
    {
        double externalWidth; //< Description
        double internalWidth; //< Description
    };

    OVDBUTIL_API trimesh::TriMesh* generateBoolcom(ovdbutil::TwoTrimesh* mesh, const int type, const BooleanParameter& param, ccglobal::Tracer* tracer=nullptr);
    OVDBUTIL_API trimesh::TriMesh* generateDilationcom(trimesh::TriMesh* mesh, const int param, ccglobal::Tracer* tracer = nullptr);
    OVDBUTIL_API trimesh::TriMesh* generateErosioncom(trimesh::TriMesh* mesh, const int param, ccglobal::Tracer* tracer = nullptr);
    OVDBUTIL_API trimesh::TriMesh* voroniaBoolcom(std::vector<trimesh::TriMesh>* mesh, const int type, const ovdbutil::BooleanParameter& param, ccglobal::Tracer* tracer=nullptr);
    OVDBUTIL_API trimesh::TriMesh* voroniaBoolcomInputLines(std::vector<edge> &edges);
}

#endif // OVDBUTIL_HOLLOWING_1650957593077_H