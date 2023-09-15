#include "ray.h"

namespace mmesh
{
	Ray::Ray()
		:start(0.0f, 0.0f, 0.0f)
		,dir(0.0f, 0.0f, 1.0f)
	{

	}

	Ray::Ray(const trimesh::vec3& vstart)
		:start(vstart)
		,dir(0.0f, 0.0f, 1.0f)
	{

	}

	Ray::Ray(const trimesh::vec3& vstart, const trimesh::vec3& ndir)
		:start(vstart)
		,dir(ndir)
	{

	}

	Ray::~Ray()
	{
	}
}
