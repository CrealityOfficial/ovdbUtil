#ifndef OVDBUTIL_HOLLOWING_1650957593077_H
#define OVDBUTIL_HOLLOWING_1650957593077_H
#include "ovdbutil/interface.h"
#include <vector>
#include "trimesh2/TriMesh.h"

/*! \file hollowing.h
    \brief A Documented file ���ͷ�ļ�������Ϊ�ӿں����ṩ���ļ�ͷ.

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

    enum PRECESION {
        COARSE =0x1,
        NORMAL=0x2,
        ELABORATE=0x3,
        EXTREME=0x4
    };
	
    /**
       *��ǵ��ڲ�������ò����ṹ����HollowingParameter�ṹ���ӽṹ��generateInfill��������õ�

       *\param enable ���ʹ�ܣ��Ƿ����
       *\param filltype ��ʱ����, ĿǰĬ��������Բ��
       *\param fillRadius ���İ뾶
       *\param fillratio ���ı���
       *\param fillLenMin ������С����
       *\param gridSizeMin MarchingCube������С��С
       *\param gridSize MarchingCube���ش�С

       *\todo ���ط���ϸ�㷨���̻�ֵ�ü����о����Ƚ϶��ѧ��֪ʶ
    */

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

  
    /**
       *��ǵĻ��������Ľṹ
       *\param min_thickness ��Ǻ��
       *\param quality �����ط����������Ľ���ľ���
       *\param closing_distance  ���������ط�interiorBandWidth�Ĵ�С
       *\param voxel_size_inout_range  ���������ط���voxel scale 
       *\param voxel_size ���������ط���voxel��С
       *\todo ���ط���ϸ�㷨���̻�ֵ�ü����о����Ƚ϶��ѧ��֪ʶ
   */
   
    struct HollowingParameter
    {
        double min_thickness = 1.0;
        double quality = 0.5;
        double closing_distance = 0.0;
        double voxel_size_inout_range = 1.0;
        double voxel_size = 1.0;
		
        int precision = NORMAL;
        bool remain_main_shell = true;
        bool filter_shell=false;
        double filter_tiny_shell = 5.0;

		INNER_FILL_CONFIG fill_config;
    };

    /// @brief  generateInterior ��������һ�����ǳ�ǵķ�������ȥ����֮����ڲ�mesh
    /// @details �������ط���levelset�ع��㷨ʵʩ
    /// 
    OVDBUTIL_API trimesh::TriMesh* generateInterior(trimesh::TriMesh* mesh, 
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    /// @brief  �պ������ڲ��޴���
    OVDBUTIL_API trimesh::TriMesh* hollowMesh(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

 /**
    *���ڽ�һ��mesh��stl���ļ����룬 �ٳ��
    *\param TriMesh ����ģ��
    *\param HollowingParameter
    *\param Tracer
    *\return TriMesh ���ģ��
    *\todo �˽ӿ��г�Ǻ����Ĺ��ܣ��Ƿ�����ڲ���䣬��INNER_FILL_CONFIG�ṹ��enable����
*/
///   @warning ������ط���ı�ԭ��ģ�;ֲ�mesh�ṹ��ʹ֮����Բ�����ط����ʾ����ģ����ӻ�����Կ���������Ƚϵ�ϸ�顣
///     ����һ������hollowMeshAndFill��ǵ�API�ӿڵĺ�����������generateInterior�����������mesh��ԭ���������mesh������ʵ���ϲ���һ���к�ȵĿ�
    OVDBUTIL_API trimesh::TriMesh* hollowMeshAndFill(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    /// @brief generateInfill������Ϊ���溯���Ļ�������������ɳ��֮���ڳ�յĿռ��������ϸС��Բ����Բ���Ĳ�����INNER_FILL_CONFIG��
	OVDBUTIL_API std::vector<trimesh::TriMesh*> generateInfill(trimesh::TriMesh* mesh, const trimesh::vec3& normal,
		const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    OVDBUTIL_API trimesh::TriMesh* hollowPrecisionMeshAndFill(trimesh::TriMesh* mesh,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    OVDBUTIL_API trimesh::TriMesh* SelectFacesHollow(trimesh::TriMesh* mesh,std::vector<int>& selectfaces,
        const HollowingParameter & = HollowingParameter(), ccglobal::Tracer* tracer = nullptr);

    OVDBUTIL_API void FindShellVolume(trimesh::TriMesh* mesh, const HollowingParameter & = HollowingParameter());
    OVDBUTIL_API bool CheckConnectChunk(trimesh::TriMesh* mesh, std::vector<std::vector<int>>& chunks,std::vector<int>& block);
}

#endif // OVDBUTIL_HOLLOWING_1650957593077_H