#ifndef _MMESH_RAY_1588946909315_H
#define _MMESH_RAY_1588946909315_H
#include "trimesh2/Vec.h"

namespace mmesh
{
	class Ray
	{
	public:
		Ray();
		Ray(const trimesh::vec3& vstart);
		Ray(const trimesh::vec3& vstart, const trimesh::vec3& ndir);

		~Ray();

		trimesh::vec3 start;
		trimesh::vec3 dir;
	};
}
#endif // _qtuser_3d_RAY_1588946909315_H
