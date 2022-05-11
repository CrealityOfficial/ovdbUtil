#ifndef OVDBUTIL_HOLLOWING_1650957593077_H
#define OVDBUTIL_HOLLOWING_1650957593077_H
#include "ovdbutil/interface.h"
#include <vector>
#include "trimesh2/Vec.h"

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


    struct HollowingParameter
    {
        double min_thickness = 1.0;
        double quality = 0.5;
        double closing_distance = 0.0;
        double voxel_size_inout_range = 1.0;
        double voxel_size = 1.0;
		
		
		INNER_FILL_CONFIG fill_config;
    };
	
	

    struct TwoTrimesh
    {
        trimesh::TriMesh *m1;
        trimesh::TriMesh *m2;
    };

    OVDBUTIL_API trimesh::TriMesh* generateBoolcom(ovdbutil::TwoTrimesh* mesh, const int type, ccglobal::Tracer* tracer=nullptr);
        OVDBUTIL_API trimesh::TriMesh* generateInterior(trimesh::TriMesh* mesh, std::vector<trimesh::vec3>* supportPoints,
        const HollowingParameter& = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    OVDBUTIL_API void hollowMesh(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    //OVDBUTIL_API std::vector<trimesh::point>* generatorSupportPoint(trimesh::TriMesh* amesh, INNER_FILL_CONFIG fillConfig);
}

#endif // OVDBUTIL_HOLLOWING_1650957593077_H