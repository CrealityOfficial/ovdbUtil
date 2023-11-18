#ifndef OVDBUTIL_GRID_1593580656409_H
#define OVDBUTIL_GRID_1593580656409_H
#include "ovdbutil/interface.h"

namespace ovdbutil
{
	struct GridBuildParam
	{
		float               exteriorBandWidth = 3.0f;
		float               interiorBandWidth = 3.0f;
		double              voxel_size = 1.0f;
	};

	class GridObjectImpl;
	class OVDBUTIL_API GridObject
	{
	public:
		GridObject();
		~GridObject();

		void setInput(DBMeshPtr mesh, const GridBuildParam& param,
			ccglobal::Tracer* tracer = nullptr);

		void dump(const std::string& fileName);
		void interateLeaf();
	protected:
		void clear();
	protected:
		GridObjectImpl* impl;
	};

	typedef std::shared_ptr<GridObject> GridObjectPtr;
}
#endif // OVDBUTIL_GRID_1593580656409_H