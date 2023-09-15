#include "trimeshutil.h"
#include "trimesh2/TriMesh.h"
#include "trimesh2/TriMesh_algo.h"
#include "trimesh2/quaternion.h"

#include <unordered_map>
#include <assert.h>

#include "ccglobal/tracer.h"
#include "ccglobal/spycc.h"
#include "ccglobal/platform.h"
#include "ccglobal/log.h"
#include <ctime>

namespace mmesh
{

	void mergeTriMesh(trimesh::TriMesh* outMesh, const std::vector<trimesh::TriMesh*>& inMeshes, bool fanzhuan)
	{
		assert(outMesh);
		size_t totalVertexSize = outMesh->vertices.size();
        size_t totalUVSize = outMesh->cornerareas.size();
		size_t totalTriangleSize = outMesh->faces.size();
        
		size_t addVertexSize = 0;
		size_t addTriangleSize = 0;
        size_t addUVSize = 0;
        
		size_t meshSize = inMeshes.size();
		for (size_t i = 0; i < meshSize; ++i)
		{
			if (inMeshes.at(i))
			{
				addVertexSize += inMeshes.at(i)->vertices.size();
				addTriangleSize += inMeshes.at(i)->faces.size();
                addUVSize += inMeshes.at(i)->cornerareas.size();
			}
		}
		totalVertexSize += addVertexSize;
		totalTriangleSize += addTriangleSize;
        totalUVSize += addUVSize;

		if (addVertexSize > 0 && addTriangleSize > 0)
		{
			outMesh->vertices.reserve(totalVertexSize);
            outMesh->cornerareas.reserve(totalUVSize);
            outMesh->faces.reserve(totalTriangleSize);

            size_t startFaceIndex = outMesh->faces.size();
            size_t startVertexIndex = outMesh->vertices.size();;
            size_t startUVIndex = outMesh->cornerareas.size();
			for (size_t i = 0; i < meshSize; ++i)
			{
				trimesh::TriMesh* mesh = inMeshes.at(i);
				if (mesh)
				{
					int vertexNum = (int)mesh->vertices.size();
					int faceNum = (int)mesh->faces.size();
                    int uvNum = (int)mesh->cornerareas.size();
					if (vertexNum > 0 && faceNum > 0)
					{
						outMesh->vertices.insert(outMesh->vertices.end(), mesh->vertices.begin(), mesh->vertices.end());
                        outMesh->cornerareas.insert(outMesh->cornerareas.end(), mesh->cornerareas.begin(), mesh->cornerareas.end());
						outMesh->faces.insert(outMesh->faces.end(), mesh->faces.begin(), mesh->faces.end());

                        size_t endFaceIndex = startFaceIndex + faceNum;
						if (startVertexIndex > 0)
						{
							for (size_t ii = startFaceIndex; ii < endFaceIndex; ++ii)
							{
								trimesh::TriMesh::Face& face = outMesh->faces.at(ii);
								for (int j = 0; j < 3; ++j)
									face[j] += startVertexIndex;

								if (fanzhuan)
								{
									int t = face[1];
									face[1] = face[2];
									face[2] = t;
								}
							}
						}

						startFaceIndex += faceNum;
						startVertexIndex += vertexNum;
                        startUVIndex += uvNum;
                        
					}
				}
			}
		}
	}

	void reverseTriMesh(trimesh::TriMesh* Mesh)
	{
		for (size_t i = 0; i < Mesh->faces.size(); i++)
		{
			int temp = Mesh->faces[i].at(1);
			Mesh->faces[i].at(1) = Mesh->faces[i].at(2);
			Mesh->faces[i].at(2) = temp;
		}
	}

	void fillTriangleSoupFaceIndex(trimesh::TriMesh* mesh)
	{
		if (!mesh || mesh->faces.size() != 0)
			return;

		size_t size = mesh->vertices.size();
		if (size % 3 || size < 3)
			return;

		size /= 3;
		mesh->faces.resize(size);
		for (size_t i = 0; i < size; ++i)
		{
			trimesh::TriMesh::Face& face = mesh->faces.at(i);
			face[0] = (int)(3 * i);
			face[1] = (int)(3 * i + 1);
			face[2] = (int)(3 * i + 2);
		}
	}

	trimesh::fxform fromQuaterian(const trimesh::quaternion& q)
	{
		float x2 = q.xp * q.xp;
		float y2 = q.yp * q.yp;
		float z2 = q.zp * q.zp;
		float xy = q.xp * q.yp;
		float xz = q.xp * q.zp;
		float yz = q.yp * q.zp;
		float wx = q.wp * q.xp;
		float wy = q.wp * q.yp;
		float wz = q.wp * q.zp;


		// This calculation would be a lot more complicated for non-unit length quaternions
		// Note: The constructor of Matrix4 expects the Matrix in column-major format like expected by
		//   OpenGL
		trimesh::fxform m = trimesh::fxform(1.0f - 2.0f * (y2 + z2), 2.0f * (xy - wz), 2.0f * (xz + wy), 0.0f,
			2.0f * (xy + wz), 1.0f - 2.0f * (x2 + z2), 2.0f * (yz - wx), 0.0f,
			2.0f * (xz - wy), 2.0f * (yz + wx), 1.0f - 2.0f * (x2 + y2), 0.0f,
			0.0f, 0.0f, 0.0f, 1.0f);

		trimesh::transpose(m);
		return m;
	}
}
