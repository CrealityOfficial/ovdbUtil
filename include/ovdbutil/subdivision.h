#ifndef OVDBUTIL_SUBDIVISION_1650957593077_H
#define OVDBUTIL_SUBDIVISION_1650957593077_H
#include "ovdbutil/interface.h"
#include <vector>
#include "trimesh2/TriMesh.h"
#include "hollowing.h"

/*! \file subdivision.h
	\brief A Documented file ϸ��ͷ�ļ�������Ϊ�ӿں����ṩ���ļ�ͷ.

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
   *���ڽ�һ��mesh��stl���ļ����룬 ��ϸ��ģ�͵���������ʹ��ø���Բ��
   *\param TriMesh ����ģ��
   *\param HollowingParameter
   *\param sub ϸ��ģ�;�ϸ���̶�ֵ
   *\param Tracer
   *\return TriMesh ���ģ��
   *\todo �˽ӿڲ�����һ����������Ľ���Բ��Ч�����������д�����
*/
   ///   @warning ������ط���ı�ԭ��ģ�;ֲ�mesh�ṹ��ʹ֮����Բ�����ط����ʾ����ģ����ӻ�����Կ���������Ƚϵ�ϸ�顣
   ///   @warning ����sub��ֵԽС���㷨��ʱԽ��
   

	OVDBUTIL_API trimesh::TriMesh* subdivison(trimesh::TriMesh* mesh,
		const SubdivisonParameter & = SubdivisonParameter(), const float sub=0.25, ccglobal::Tracer* tracer = nullptr);

/**
   *���ڽ�һ��mesh��stl���ļ����룬���Ƴ�ǣ�������ȫ�ǣ�Ŀǰ�ú�������ģ��voronoi�οյ���Ŀ�У���һ���Ƚ�����Ľӿڡ�
   *\param TriMesh ����ģ��
   *\param xformf ���ط������س߶Ȳ�����ϸ��ģ�;�ϸ���̶�ֵ��ԽС��Խ��ϸ������ʱ��Խ������ǵĿǸ��ӵĺ񣬴�������Ч��
   *\param w ����ֵ�����Կǵĺ�����ں����������
   *\param Tracer
   *\return TriMesh ���ģ��
   *\todo �˽ӿ���
*/
///   @warning voronoi�οյ���Ŀ�е�ר�нӿڡ�
/// 
	OVDBUTIL_API trimesh::TriMesh* remesh(trimesh::TriMesh* mesh, const float xformf, const int w, ccglobal::Tracer* tracer);
}

#endif // OVDBUTIL_SUBDIVISION_1650957593077_H