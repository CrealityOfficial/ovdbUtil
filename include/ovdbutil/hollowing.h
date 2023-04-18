#ifndef OVDBUTIL_HOLLOWING_1650957593077_H
#define OVDBUTIL_HOLLOWING_1650957593077_H
#include "ovdbutil/interface.h"
#include <vector>
#include "trimesh2/TriMesh.h"

/*! \file hollowing.h
    \brief A Documented file 抽壳.

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
		bool enable = false;//填充使能
		int filltype = 0;
		float fillRadius = 1.0;//填充的半径
		float fillratio = 0.5;//填充的比率
		float fillLenMin = 1.0;//填充的最小长度
		float gridSizeMin = 5.0;//MarchingCube体素最小大小
		float gridSize = 5.0;//MarchingCube体素大小
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
    *用于将一个mesh的stl的文件输入， 再抽壳
    *\li 测试中。。。。。。1 测试修改了之后会不会覆盖原来的东西
    *\li 测试中。。。。。。2 测试修改了之后会不会覆盖原来的东西
    *\li 测试中。。。。。。3

    *\param TriMesh
    *\param HollowingParameter
    *\param Tracer

    *\return TriMesh

    *\todo 在xxx的此接口有抽壳和填充的功能
*/

/// @brief Implementation ImplementationImplementationImplementationImplementationImplementationImplementation
/// @details 
/// 
///     创建一个测试hollowMeshAndFill抽壳的API接口的函数
    OVDBUTIL_API trimesh::TriMesh* hollowMeshAndFill(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

	OVDBUTIL_API std::vector<trimesh::TriMesh*> generateInfill(trimesh::TriMesh* mesh, const trimesh::vec3& normal,
		const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);
}

#endif // OVDBUTIL_HOLLOWING_1650957593077_H