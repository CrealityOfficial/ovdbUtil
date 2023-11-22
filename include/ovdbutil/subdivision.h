#ifndef OVDBUTIL_SUBDIVISION_1650957593077_H
#define OVDBUTIL_SUBDIVISION_1650957593077_H
#include "ovdbutil/interface.h"
#include <vector>
#include "trimesh2/TriMesh.h"
#include "hollowing.h"

/*! \file subdivision.h
	\brief A Documented file 细分头文件，单独为接口函数提供的文件头.

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


    struct SubdivisonParameter
    {
        double min_thickness = 1.0;
        double quality = 0.5;
        double closing_distance = 0.0;
        double voxel_size_inout_range = 1.0;
        double voxel_size = 1.0;
		
		
		INNER_FILL_CONFIG fill_config;
    };

/**
   *用于将一个mesh的stl的文件输入， 再细分模型的三角网格，使变得更加圆润。
   *\param TriMesh 输入模型
   *\param HollowingParameter
   *\param sub 细分模型精细化程度值
   *\param Tracer
   *\return TriMesh 输出模型
   *\todo 此接口不像抽壳一样对立方体的角有圆润效果。具体结果有待测试
*/
   ///   @warning 这个体素法会改变原来模型局部mesh结构，使之更加圆润，体素法本质决定的，可视化后可以看到三角面比较的细碎。
   ///   @warning 参数sub数值越小，算法耗时越多
   

	OVDBUTIL_API trimesh::TriMesh* subdivison(trimesh::TriMesh* mesh,
		const SubdivisonParameter & = SubdivisonParameter(), const float sub=0.25, ccglobal::Tracer* tracer = nullptr);

/**
   *用于将一个mesh的stl的文件输入，类似抽壳，但不完全是，目前该函数用在模型voronoi镂空的项目中，是一个比较特殊的接口。
   *\param TriMesh 输入模型
   *\param xformf 体素法的体素尺度参数，细分模型精细化程度值，越小，越精细，计算时间越长。抽壳的壳更加的厚，带有膨胀效果
   *\param w 膨胀值，即对壳的厚度向内和外进行膨胀
   *\param Tracer
   *\return TriMesh 输出模型
   *\todo 此接口有
*/
///   @warning voronoi镂空的项目中的专有接口。
/// 
	OVDBUTIL_API trimesh::TriMesh* remesh(trimesh::TriMesh* mesh, const float xformf, const int w, ccglobal::Tracer* tracer);
}

#endif // OVDBUTIL_SUBDIVISION_1650957593077_H