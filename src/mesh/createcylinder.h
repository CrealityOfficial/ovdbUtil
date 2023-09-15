#ifndef MMESH_CREATECYLINDER_1619866821989_H
#define MMESH_CREATECYLINDER_1619866821989_H
#include "trimesh2/TriMesh.h"

namespace mmesh
{
	trimesh::TriMesh* createSoupCylinder(int count, float _radius, float _height, bool offsetOnZero = false);
	trimesh::TriMesh* createSoupCylinder(int count, float _radius, float _height,
		const trimesh::vec3& centerPoint, const trimesh::vec3& normal);
}

#endif // MMESH_CREATECYLINDER_1619866821989_H