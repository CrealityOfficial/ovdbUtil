#include "ovdbutil/grid.h"

#include "util/gridhelper.h"
#include "util/tracer.h"

namespace ovdbutil
{
	class GridObjectImpl
	{
	public:
		GridObjectImpl() {}
		~GridObjectImpl() {}

		DBMeshPtr mesh;
		FloatGridPtr grid;
	};

	GridObject::GridObject()
		:impl(new GridObjectImpl())
	{

	}

	GridObject::~GridObject()
	{

	}

	void GridObject::setInput(DBMeshPtr mesh, const GridBuildParam& param,
		ccglobal::Tracer* tracer)
	{
		clear();

		impl->mesh = mesh;
		TracerInterrupter interrupter(tracer);

		openvdb::math::Transform tr;
		impl->grid = mesh_to_grid(interrupter, mesh.get(), tr,
			param.exteriorBandWidth, param.interiorBandWidth, param.voxel_size);
	}

	void GridObject::dump(const std::string& fileName)
	{
	}

	void GridObject::interateLeaf()
	{
		if (!impl->grid)
			return;
	}

	void GridObject::clear()
	{
		impl->mesh = nullptr;
		impl->grid = nullptr;
	}
}