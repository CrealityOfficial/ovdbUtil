#ifndef OVDBUTIL_SUBDIVISION_1650957593077_H
#define OVDBUTIL_SUBDIVISION_1650957593077_H
#include "ovdbutil/interface.h"
#include <vector>
#include "trimesh2/TriMesh.h"

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
	typedef struct INNER_FILL_CONFIG
	{
		bool enable = false;//���ʹ��
		int filltype = 0;
		float fillRadius = 1.0;//���İ뾶
		float fillratio = 0.5;//���ı���
		float fillLenMin = 1.0;//������С����
		float gridSizeMin = 5.0;//MarchingCube������С��С
		float gridSize = 5.0;//MarchingCube���ش�С
	}sINNER_FILL_CFG;


    struct SubdivisonParameter
    {
        double min_thickness = 1.0;
        double quality = 0.5;
        double closing_distance = 0.0;
        double voxel_size_inout_range = 1.0;
        double voxel_size = 1.0;
		
		
		INNER_FILL_CONFIG fill_config;
    };



	OVDBUTIL_API trimesh::TriMesh* subdivison(trimesh::TriMesh* mesh,
		const SubdivisonParameter & = SubdivisonParameter(), const float sub=0.25, ccglobal::Tracer* tracer = nullptr);
	OVDBUTIL_API trimesh::TriMesh* remesh(trimesh::TriMesh* mesh, const float xformf, const int w, ccglobal::Tracer* tracer);
}

#endif // OVDBUTIL_SUBDIVISION_1650957593077_H