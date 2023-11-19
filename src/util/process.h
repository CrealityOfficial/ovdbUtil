#ifndef PROCESS_H_H
#define PROCESS_H_H
#include <openvdb/tree/LeafManager.h>

namespace ovdbutil
{
    namespace util {

        /// Helper class used internally by processTypedGrid()
        template<typename GridType, typename OpType, bool IsConst/*=false*/>
        struct GridProcessor {
            static inline void call(OpType& op, openvdb::GridBase::Ptr grid) {
                op.template operator() < GridType > (openvdb::gridPtrCast<GridType>(grid));
            }
        };

        /// Helper class used internally by processTypedGrid()
        template<typename GridType, typename OpType>
        struct GridProcessor<GridType, OpType, /*IsConst=*/true> {
            static inline void call(OpType& op, openvdb::GridBase::ConstPtr grid) {
                op.template operator() < GridType > (openvdb::gridConstPtrCast<GridType>(grid));
            }
        };


        /// Helper function used internally by processTypedGrid()
        template<typename GridType, typename OpType, typename GridPtrType>
        inline void
            doProcessTypedGrid(GridPtrType grid, OpType& op)
        {
            GridProcessor<GridType, OpType,
                std::is_const<typename GridPtrType::element_type>::value>::call(op, grid);
        }


        ////////////////////////////////////////


        /// @brief Utility function that, given a generic grid pointer,
        /// calls a functor on the fully-resolved grid
        ///
        /// Usage:
        /// @code
        /// struct PruneOp {
        ///     template<typename GridT>
        ///     void operator()(typename GridT::Ptr grid) const { grid->tree()->prune(); }
        /// };
        ///
        /// processTypedGrid(myGridPtr, PruneOp());
        /// @endcode
        ///
        /// @return @c false if the grid type is unknown or unhandled.
        template<typename GridPtrType, typename OpType>
        bool
            processTypedGrid(GridPtrType grid, OpType& op)
        {
            using namespace openvdb;
            if (grid->template isType<BoolGrid>())        doProcessTypedGrid<BoolGrid>(grid, op);
            else if (grid->template isType<FloatGrid>())  doProcessTypedGrid<FloatGrid>(grid, op);
            else if (grid->template isType<DoubleGrid>()) doProcessTypedGrid<DoubleGrid>(grid, op);
            else if (grid->template isType<Int32Grid>())  doProcessTypedGrid<Int32Grid>(grid, op);
            else if (grid->template isType<Int64Grid>())  doProcessTypedGrid<Int64Grid>(grid, op);
            else if (grid->template isType<Vec3IGrid>())  doProcessTypedGrid<Vec3IGrid>(grid, op);
            else if (grid->template isType<Vec3SGrid>())  doProcessTypedGrid<Vec3SGrid>(grid, op);
            else if (grid->template isType<Vec3DGrid>())  doProcessTypedGrid<Vec3DGrid>(grid, op);
            else if (grid->template isType<points::PointDataGrid>()) {
                doProcessTypedGrid<points::PointDataGrid>(grid, op);
            }
            else return false;
            return true;
        }


        /// @brief Utility function that, given a generic grid pointer, calls
        /// a functor on the fully-resolved grid, provided that the grid's
        /// voxel values are scalars
        template<typename GridPtrType, typename OpType>
        bool
            processTypedScalarGrid(GridPtrType grid, OpType& op)
        {
            using namespace openvdb;
            if (grid->template isType<FloatGrid>())       doProcessTypedGrid<FloatGrid>(grid, op);
            else if (grid->template isType<DoubleGrid>()) doProcessTypedGrid<DoubleGrid>(grid, op);
            else if (grid->template isType<Int32Grid>())  doProcessTypedGrid<Int32Grid>(grid, op);
            else if (grid->template isType<Int64Grid>())  doProcessTypedGrid<Int64Grid>(grid, op);
            else return false;
            return true;
        }


        /// @brief Utility function that, given a generic grid pointer, calls
        /// a functor on the fully-resolved grid, provided that the grid's
        /// voxel values are scalars or PointIndex objects
        template<typename GridPtrType, typename OpType>
        bool
            processTypedScalarOrPointDataGrid(GridPtrType grid, OpType& op)
        {
            using namespace openvdb;
            if (processTypedScalarGrid(grid, op)) return true;
            if (grid->template isType<points::PointDataGrid>()) {
                doProcessTypedGrid<points::PointDataGrid>(grid, op);
                return true;
            }
            return false;
        }


        /// @brief Utility function that, given a generic grid pointer, calls
        /// a functor on the fully-resolved grid, provided that the grid's
        /// voxel values are vectors
        template<typename GridPtrType, typename OpType>
        bool
            processTypedVectorGrid(GridPtrType grid, OpType& op)
        {
            using namespace openvdb;
            if (grid->template isType<Vec3IGrid>())       doProcessTypedGrid<Vec3IGrid>(grid, op);
            else if (grid->template isType<Vec3SGrid>())  doProcessTypedGrid<Vec3SGrid>(grid, op);
            else if (grid->template isType<Vec3DGrid>())  doProcessTypedGrid<Vec3DGrid>(grid, op);
            else return false;
            return true;
        }

        template<class TreeType>
        class MinMaxVoxel
        {
        public:
            using LeafArray = openvdb::tree::LeafManager<TreeType>;
            using ValueType = typename TreeType::ValueType;

            // LeafArray = openvdb::tree::LeafManager<TreeType> leafs(myTree)
            MinMaxVoxel(LeafArray&);

            void runParallel();
            void runSerial();

            const ValueType& minVoxel() const { return mMin; }
            const ValueType& maxVoxel() const { return mMax; }

            inline MinMaxVoxel(const MinMaxVoxel<TreeType>&, tbb::split);
            inline void operator()(const tbb::blocked_range<size_t>&);
            inline void join(const MinMaxVoxel<TreeType>&);

        private:
            LeafArray& mLeafArray;
            ValueType mMin, mMax;
        };


        template <class TreeType>
        MinMaxVoxel<TreeType>::MinMaxVoxel(LeafArray& leafs)
            : mLeafArray(leafs)
            , mMin(std::numeric_limits<ValueType>::max())
            , mMax(std::numeric_limits<ValueType>::lowest())
        {
        }


        template <class TreeType>
        inline
            MinMaxVoxel<TreeType>::MinMaxVoxel(const MinMaxVoxel<TreeType>& rhs, tbb::split)
            : mLeafArray(rhs.mLeafArray)
            , mMin(std::numeric_limits<ValueType>::max())
            , mMax(std::numeric_limits<ValueType>::lowest())
        {
        }


        template <class TreeType>
        void
            MinMaxVoxel<TreeType>::runParallel()
        {
            tbb::parallel_reduce(mLeafArray.getRange(), *this);
        }


        template <class TreeType>
        void
            MinMaxVoxel<TreeType>::runSerial()
        {
            (*this)(mLeafArray.getRange());
        }


        template <class TreeType>
        inline void
            MinMaxVoxel<TreeType>::operator()(const tbb::blocked_range<size_t>& range)
        {
            typename TreeType::LeafNodeType::ValueOnCIter iter;

            for (size_t n = range.begin(); n < range.end(); ++n) {
                iter = mLeafArray.leaf(n).cbeginValueOn();
                for (; iter; ++iter) {
                    const ValueType value = iter.getValue();
                    mMin = std::min(mMin, value);
                    mMax = std::max(mMax, value);
                }
            }
        }


        template <class TreeType>
        inline void
            MinMaxVoxel<TreeType>::join(const MinMaxVoxel<TreeType>& rhs)
        {
            mMin = std::min(mMin, rhs.mMin);
            mMax = std::max(mMax, rhs.mMax);
        }

    } // namespace util
}

#endif // PROCESS_H_H
