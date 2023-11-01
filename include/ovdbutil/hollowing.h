#ifndef OVDBUTIL_HOLLOWING_1650957593077_H
#define OVDBUTIL_HOLLOWING_1650957593077_H
#include "ovdbutil/interface.h"
#include <vector>
#include "trimesh2/TriMesh.h"


/*! \file hollowing.h
    \brief A Documented file 抽壳头文件，单独为接口函数提供的文件头.

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
    enum PRECISION
    {
        EXTREMITY = 4,
        EXQUISITE = 3,
        NORMAL = 2,
        CRUDE = 1,
    };
    /**
       *抽壳的内部填充配置参数结构，是HollowingParameter结构的子结构，generateInfill这个函数用到

       *\param enable 填充使能，是否填充
       *\param filltype 暂时无用, 目前默认填充的是圆柱
       *\param fillRadius 填充的半径
       *\param fillratio 填充的比率
       *\param fillLenMin 填充的最小长度
       *\param gridSizeMin MarchingCube体素最小大小
       *\param gridSize MarchingCube体素大小

       *\todo 体素法详细算法过程还值得继续研究，比较多的学术知识
    */

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

  
    /**
       *抽壳的基本参数的结构
       *\param min_thickness 抽壳厚度
       *\param quality 用体素法产生的最后的结果的精度
       *\param closing_distance  决定了体素法interiorBandWidth的大小
       *\param voxel_size_inout_range  决定了体素法的voxel scale 
       *\param voxel_size 决定了体素法的voxel大小
       *\todo 体素法详细算法过程还值得继续研究，比较多的学术知识
   */
   
    struct HollowingParameter
    {
        double min_thickness = 1.0;
        double quality = 0.5;
        double closing_distance = 0.0;
        double voxel_size_inout_range = 1.0;
        double voxel_size = 1.0;
		
        int precision = CRUDE;
		
		INNER_FILL_CONFIG fill_config;
    };

   

    /// @brief  generateInterior 像函数名字一样，是抽壳的反操作，去掉壳之后的内部mesh
    /// @details 基于体素法的levelset重构算法实施
    /// 
    OVDBUTIL_API trimesh::TriMesh* generateInterior(trimesh::TriMesh* mesh, 
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    /// @brief  启用这个函数当作控制精细程度的接口函数
    OVDBUTIL_API trimesh::TriMesh* hollowMesh(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);
    
 /**
    *用于将一个mesh的stl的文件输入， 再抽壳
    *\param TriMesh 输入模型
    *\param HollowingParameter
    *\param Tracer
    *\return TriMesh 输出模型
    *\todo 此接口有抽壳和填充的功能，是否产生内部填充，由INNER_FILL_CONFIG结构的enable控制
*/
///   @warning 这个体素法会改变原来模型局部mesh结构，使之更加圆润，体素法本质决定的，可视化后可以看到三角面比较的细碎。
///     创建一个测试hollowMeshAndFill抽壳的API接口的函数，真正的generateInterior函数抽壳里面mesh和原来的输入的mesh，这样实际上才是一个有厚度的壳
    OVDBUTIL_API trimesh::TriMesh* hollowMeshAndFill(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    /// @brief generateInfill函数做为上面函数的基础函数，在完成抽壳之后，在抽空的空间区域，填充细小的圆柱，圆柱的参数在INNER_FILL_CONFIG中
	OVDBUTIL_API std::vector<trimesh::TriMesh*> generateInfill(trimesh::TriMesh* mesh, const trimesh::vec3& normal,
		const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    OVDBUTIL_API trimesh::TriMesh* hollowPrecisionMeshAndFill(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);
}

#endif // OVDBUTIL_HOLLOWING_1650957593077_H