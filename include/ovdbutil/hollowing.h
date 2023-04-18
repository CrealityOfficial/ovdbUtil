#ifndef OVDBUTIL_HOLLOWING_1650957593077_H
#define OVDBUTIL_HOLLOWING_1650957593077_H
#include "ovdbutil/interface.h"
#include <vector>
#include "trimesh2/TriMesh.h"

/*! \file hollowing.h
    \brief A Documented file ���.

    Details.
*/

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
	
	/// @brief INNER_FILL_CONFIG

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

  
    /// @brief HollowingParameter
   
    struct HollowingParameter
    {
        double min_thickness = 1.0;
        double quality = 0.5;
        double closing_distance = 0.0;
        double voxel_size_inout_range = 1.0;
        double voxel_size = 1.0;
		
		
		INNER_FILL_CONFIG fill_config;
    };

    /// @brief  generateInterior Implementation ImplementationImplementationImplementationImplementationImplementationImplementation
    /// @details 
    /// 
    OVDBUTIL_API trimesh::TriMesh* generateInterior(trimesh::TriMesh* mesh, std::vector<trimesh::vec3>* supportPoints,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);
    /// @brief  hollowMesh Implementation ImplementationImplementationImplementationImplementationImplementationImplementation
    /// @details 
    /// 
    OVDBUTIL_API void hollowMesh(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

 /**
    *���ڽ�һ��mesh��stl���ļ����룬 �ٳ��
    *\li �����С�����������1 �����޸���֮��᲻�Ḳ��ԭ���Ķ���
    *\li �����С�����������2 �����޸���֮��᲻�Ḳ��ԭ���Ķ���
    *\li �����С�����������3

    *\param TriMesh
    *\param HollowingParameter
    *\param Tracer

    *\return TriMesh

    *\todo ��xxx�Ĵ˽ӿ��г�Ǻ����Ĺ���
*/

/// @brief Implementation ImplementationImplementationImplementationImplementationImplementationImplementation
/// @details 
/// 
///     ����һ������hollowMeshAndFill��ǵ�API�ӿڵĺ���
    OVDBUTIL_API trimesh::TriMesh* hollowMeshAndFill(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

	OVDBUTIL_API std::vector<trimesh::TriMesh*> generateInfill(trimesh::TriMesh* mesh, const trimesh::vec3& normal,
		const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);
}

#endif // OVDBUTIL_HOLLOWING_1650957593077_H