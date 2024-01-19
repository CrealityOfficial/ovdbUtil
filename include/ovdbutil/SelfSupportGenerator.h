#pragma once

#include <openvdb/math/Vec3.h>
#include <openvdb/math/Coord.h>
#include <openvdb/tools/ValueTransformer.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/GridTransformer.h>
#include "VOXELTOOL/VoxelBody.hpp"

namespace ovdbutil
{
    class SelfSupportGenerator
    {

        using ValueT = typename openvdb::FloatGrid::ValueType;

    public:
        SelfSupportGenerator(openvdb::FloatGrid::Ptr& gridptr, const double& offset) : gridptr(gridptr), offset(offset)
        {
            openvdb::Coord boxStart = gridptr->evalActiveVoxelBoundingBox().getStart();
            openvdb::Coord boxEnd = gridptr->evalActiveVoxelBoundingBox().getEnd();
            xStart = boxStart.x();
            yStart = boxStart.y();
            zStart = boxStart.z();
            xEnd = boxEnd.x();
            yEnd = boxEnd.y();
            zEnd = boxEnd.z();
            background = gridptr->background();
        }

        void toVoxelEnumType();
        void setBackGround(const int& x, const int& y, const int& z, const VOXELTOOL::VoxelBody& voxelBody);
        void generate();
        openvdb::Coord toCoord(const int& x, const int& y, const int& z)
        {
            return openvdb::Coord(x + xStart, y + yStart, zEnd - z);
        }
        void setValue(const int& x, const int& y, const int& z, const ValueT& v)
        {
            typename openvdb::FloatGrid::Accessor accessor = gridptr->getAccessor();
            openvdb::Coord coord = toCoord(x, y, z);
            accessor.setValue(coord, v);
        }

    private:
        openvdb::FloatGrid::Ptr& gridptr;                                          // openvdb的指针
        std::vector<std::vector<std::vector<VOXELTOOL::VoxelEnumType>>> voxelData; // openvdb转换成的体素数据
        int xStart;
        int xEnd;
        int yStart;
        int yEnd;
        int zStart;
        int zEnd;
        double offset;
        ValueT background;                                                          // 背景值存储
    };

}