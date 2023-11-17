#ifndef OVDBUTIL_GRID_1593580656409_H
#define OVDBUTIL_GRID_1593580656409_H
#include "ovdbutil/interface.h"

namespace ovdbutil
{
	class GridObjectImpl;
	class OVDBUTIL_API GridObject
	{
	public:
		GridObject();
		~GridObject();

		void setInput(DBMeshPtr mesh);
	protected:
		void clear();
	protected:
		GridObjectImpl* impl;
	};
}
#endif // OVDBUTIL_GRID_1593580656409_H