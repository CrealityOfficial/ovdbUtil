#include "createcylinder.h"
#include "trimesh2/TriMesh_algo.h"
#include "trimesh2/quaternion.h"
#include "trimeshutil.h"

#define  PI 3.141592 

namespace mmesh
{
	trimesh::TriMesh* createSoupCylinder(int count, float _radius, float _height, bool offsetOnZero)
	{
		int  trianglesCount = count;
		if (trianglesCount < 3)
			trianglesCount = 3;
		trimesh::TriMesh* cylinderMesh = new trimesh::TriMesh();
		float radian = 2.0f * PI / (float)trianglesCount;//4¶È»»Ëã³É»¡¶È

		trimesh::vec3 top(0.0f, 0.0f, _height / 2.0f);
		trimesh::vec3 bottom(0.0f, 0.0f, - _height / 2.0f);
		trimesh::vec3 offset(0.0f, 0.0f, _height);
		std::vector<trimesh::vec3> circles(trianglesCount, trimesh::vec3());
		for (int i = 0; i < trianglesCount; ++i)
			circles[i] = trimesh::vec3(_radius * cos(i * radian),
				_radius * sin(i * radian), -_height / 2.0f);
		circles.push_back(circles[0]);
		cylinderMesh->vertices.resize(12 * trianglesCount);
		int index = 0;

		auto f = [&index, &cylinderMesh](const trimesh::vec3& v1, const trimesh::vec3& v2, const trimesh::vec3& v3) {
			cylinderMesh->vertices[index++] = v1;
			cylinderMesh->vertices[index++] = v2;
			cylinderMesh->vertices[index++] = v3;
		};
		for (int i = 0; i < trianglesCount; ++i)
		{
			f(top, circles[i] + offset, circles[i + 1] + offset);
			f(circles[i] + offset, circles[i], circles[i + 1] + offset);
			f(circles[i + 1] + offset, circles[i], circles[i + 1]);
			f(bottom, circles[i + 1], circles[i]);
		}

		fillTriangleSoupFaceIndex(cylinderMesh);

		if(offsetOnZero)
			trimesh::apply_xform(cylinderMesh, trimesh::xform::trans(0.0f, 0.0f, _height / 2.0f));
		return cylinderMesh;
	}

	trimesh::TriMesh* createSoupCylinder(int count, float _radius, float _height,
		const trimesh::vec3& centerPoint, const trimesh::vec3& normal)
	{
		trimesh::TriMesh* mesh = createSoupCylinder(count, _radius, _height);

		const trimesh::vec3 cyOriginNormal(0.0f, 0.0f, 1.0f);
		trimesh::quaternion q = trimesh::quaternion::rotationTo(normal, cyOriginNormal);
		trimesh::fxform xf = trimesh::fxform::trans(centerPoint) * mmesh::fromQuaterian(q);
		trimesh::apply_xform(mesh, trimesh::xform(xf));

		return mesh;
	}
}