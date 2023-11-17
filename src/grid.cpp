#include "ovdbutil/grid.h"

namespace ovdbutil
{
	class GridObjectImpl
	{
	public:
		GridObjectImpl() {}
		~GridObjectImpl() {}

		DBMeshPtr mesh;
	};

	GridObject::GridObject()
		:impl(new GridObjectImpl())
	{

	}

	GridObject::~GridObject()
	{

	}

	void GridObject::setInput(DBMeshPtr mesh)
	{

	}

	void GridObject::clear()
	{

	}
}