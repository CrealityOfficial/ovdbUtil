#ifndef MMESH_TRIMESHUTIL_1602590289222_H
#define MMESH_TRIMESHUTIL_1602590289222_H
#include <vector>
#include <fstream>
#include <memory>

#include "trimesh2/TriMesh.h"
#include "trimesh2/XForm.h"
#include "trimesh2/Box.h"
#include "trimesh2/quaternion.h"

typedef std::shared_ptr<trimesh::TriMesh> TriMeshPointer;

namespace ccglobal
{
	class Tracer;
}

namespace mmesh
{
	void mergeTriMesh(trimesh::TriMesh* outMesh, const std::vector<trimesh::TriMesh*>& inMeshes, bool fanzhuan = false);
	void reverseTriMesh(trimesh::TriMesh* Mesh);

	void fillTriangleSoupFaceIndex(trimesh::TriMesh* mesh);
	trimesh::fxform fromQuaterian(const trimesh::quaternion& q);
}

#endif // MMESH_TRIMESHUTIL_1602590289222_H