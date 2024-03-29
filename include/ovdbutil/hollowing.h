#ifndef OVDBUTIL_HOLLOWING_1650957593077_H
#define OVDBUTIL_HOLLOWING_1650957593077_H
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
		bool enable = false;//填充使能
		int filltype = 0;
		float fillRadius = 1.0;//填充的半径
		float fillratio = 0.5;//填充的比率
		float fillLenMin = 1.0;//填充的最小长度
		float gridSizeMin = 5.0;//MarchingCube体素最小大小
		float gridSize = 5.0;//MarchingCube体素大小
	}sINNER_FILL_CFG;


    struct HollowingParameter
    {
        double min_thickness = 1.0;
        double quality = 0.5;
        double closing_distance = 0.0;
        double voxel_size_inout_range = 1.0;
        double voxel_size = 1.0;
		
		
		INNER_FILL_CONFIG fill_config;
    };

        
    OVDBUTIL_API trimesh::TriMesh* generateInterior(trimesh::TriMesh* mesh, std::vector<trimesh::vec3>* supportPoints,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);
    OVDBUTIL_API void hollowMesh(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    OVDBUTIL_API trimesh::TriMesh* hollowMeshAndFill(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

	OVDBUTIL_API std::vector<trimesh::TriMesh*> generateInfill(trimesh::TriMesh* mesh, const trimesh::vec3& normal,
		const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);
}

#endif // OVDBUTIL_HOLLOWING_1650957593077_H