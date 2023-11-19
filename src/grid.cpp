#include "ovdbutil/grid.h"

#include "util/gridhelper.h"
#include "util/tracer.h"
#include "util/process.h"

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
		tr.preScale((double)param.scale);
		impl->grid = mesh_to_grid(interrupter, mesh.get(), tr,
			param.exteriorBandWidth, param.interiorBandWidth, param.voxel_size);
	}

	void GridObject::traitBoxes(VDBTopoDetail& detail)
	{
		detail.L1Boxes.clear();
		detail.L2Boxes.clear();
		detail.LeafBoxes.clear();

		if (impl->grid)
		{
			FloatGridPtr grid = impl->grid;
			typedef openvdb::FloatGrid::ValueType ValueType;
			typedef openvdb::FloatGrid::TreeType TreeType;
			typedef openvdb::FloatGrid::TreeType::NodeCIter NodeCIter;

			const bool isLevelSetGrid = grid->getGridClass() == openvdb::GRID_LEVEL_SET;
			TreeType& tree = grid->tree();

			detail.voxel_size = grid->voxelSize();
			detail.leaf_count = tree.leafCount();
			detail.no_leaf_node_count = tree.nonLeafCount();

			openvdb::Index64 nodeCount = (openvdb::Index64)tree.leafCount()
				+ (openvdb::Index64)tree.nonLeafCount();
			const openvdb::Index64 N = nodeCount * 8 * 3;

			openvdb::CoordBBox bbox;
			for (NodeCIter iter = grid->tree().cbeginNode(); iter; ++iter)
			{
				iter.getBoundingBox(bbox);

				// Nodes are rendered as cell-centered
				const openvdb::Vec3d vmin(bbox.min().x() - 0.5, bbox.min().y() - 0.5, bbox.min().z() - 0.5);
				const openvdb::Vec3d vmax(bbox.max().x() + 0.5, bbox.max().y() + 0.5, bbox.max().z() + 0.5);
				
				trimesh::box3 b;
				b += to_vec3f(grid->indexToWorld(vmin));
				b += to_vec3f(grid->indexToWorld(vmax));

				const int level = iter.getLevel();
				if (level == 0)
					detail.LeafBoxes.push_back(b);
				else if (level == 1)
					detail.L2Boxes.push_back(b);
				else if (level == 2)
					detail.L1Boxes.push_back(b);
				else
					detail.RootBoxes.push_back(b);
			}

			openvdb::tree::LeafManager<const TreeType> leafs(tree);

			{
				util::MinMaxVoxel<const TreeType> minmax(leafs);
				minmax.runParallel();
				detail.minValue = minmax.minVoxel();
				detail.maxValue = minmax.maxVoxel();
			}

            openvdb::Index64 voxelsPerLeaf = TreeType::LeafNodeType::NUM_VOXELS;

            // Level set rendering
            if (tree.activeLeafVoxelCount() > 26000000) {
                voxelsPerLeaf = std::max<openvdb::Index64>(1, (26000000 / tree.leafCount()));
            }

            std::vector<size_t> indexMap(leafs.leafCount());
            size_t voxelCount = 0;
            for (openvdb::Index64 l = 0, L = leafs.leafCount(); l < L; ++l) {
                indexMap[l] = voxelCount;
                voxelCount += std::min(leafs.leaf(l).onVoxelCount(), voxelsPerLeaf);
            }

			detail.voxel_count = voxelCount;

			using ValueOnCIter = typename TreeType::LeafNodeType::ValueOnCIter;

			for (size_t n = 0; n < leafs.leafCount(); ++n) {
				ValueOnCIter it = leafs.leaf(n).cbeginValueOn();

				openvdb::Index64 activeVoxels = leafs.leaf(n).onVoxelCount();

				if (activeVoxels <= voxelsPerLeaf) {

					for (; it; ++it) {
						detail.voxels.push_back(to_vec3f(grid->indexToWorld(it.getCoord())));
					}

				}
				else if (1 == voxelsPerLeaf) {

					detail.voxels.push_back(to_vec3f(grid->indexToWorld(it.getCoord())));
				}
				else {

					std::vector<openvdb::Coord> coords;
					coords.reserve(static_cast<size_t>(activeVoxels));
					for (; it; ++it) { coords.push_back(it.getCoord()); }

					detail.voxels.push_back(to_vec3f(grid->indexToWorld(coords[0])));
					detail.voxels.push_back(to_vec3f(grid->indexToWorld(coords[static_cast<size_t>(activeVoxels - 1)])));

					openvdb::Index64 r = openvdb::Index64(std::floor(double(voxelsPerLeaf) / double(activeVoxels)));
					for (openvdb::Index64 i = 1, I = voxelsPerLeaf - 2; i < I; ++i) {
						detail.voxels.push_back(to_vec3f(grid->indexToWorld(coords[static_cast<size_t>(i * r)])));
					}
				}
			}
		}
	}

	void GridObject::clear()
	{
		impl->mesh = nullptr;
		impl->grid = nullptr;
	}
}