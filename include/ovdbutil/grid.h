#ifndef OVDBUTIL_GRID_1593580656409_H
#define OVDBUTIL_GRID_1593580656409_H
#include "ovdbutil/interface.h"

namespace ovdbutil
{
	struct VDBTopoDetail
	{
		std::vector<trimesh::box3> RootBoxes;
		std::vector<trimesh::box3> L1Boxes;
		std::vector<trimesh::box3> L2Boxes;
		std::vector<trimesh::box3> LeafBoxes;

		float minValue = -1.0f;
		float maxValue = 1.0f;
		std::vector<trimesh::vec3> voxels;

		trimesh::vec3 voxel_size;
		int leaf_count = 0;
		int no_leaf_node_count = 0;
		int voxel_count = 0;
	};

	struct GridBuildParam
	{
		float               exteriorBandWidth = 3.0f;
		float               interiorBandWidth = 3.0f;
		double              voxel_size = 1.0f;
		float				scale = 1.0f;
	};

	class GridObjectImpl;
	class OVDBUTIL_API GridObject
	{
	public:
		GridObject();
		~GridObject();

		void setInput(DBMeshPtr mesh, const GridBuildParam& param,
			ccglobal::Tracer* tracer = nullptr);

		void traitBoxes(VDBTopoDetail& detail);
	protected:
		void clear();
	protected:
		GridObjectImpl* impl;
	};

	typedef std::shared_ptr<GridObject> GridObjectPtr;
}
#endif // OVDBUTIL_GRID_1593580656409_H