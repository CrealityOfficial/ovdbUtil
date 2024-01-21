#include "ovdbutil/SelfSupportGenerator.h"

namespace ovdbutil
{
    void SelfSupportGenerator::toVoxelEnumType()
    {
        typename openvdb::FloatGrid::Accessor accessor = gridptr->getAccessor();
        openvdb::Coord boxStart = gridptr->evalActiveVoxelBoundingBox().getStart();
        openvdb::Coord boxEnd = gridptr->evalActiveVoxelBoundingBox().getEnd();
        voxelData.reserve(boxEnd.z() - boxStart.z() + 1);
        for (int z = boxEnd.z(); z >= boxStart.z(); z--)
        {
            std::vector<std::vector<VOXELTOOL::VoxelEnumType>> zData;
            zData.reserve(boxEnd.z() - boxStart.z() + 1);
            for (int y = boxStart.y(); y <= boxEnd.y(); ++y)
            {
                std::vector<VOXELTOOL::VoxelEnumType> yData;
                yData.reserve(boxEnd.x() - boxStart.x() + 1);
                for (int x = boxStart.x(); x <= boxEnd.x(); ++x)
                {
                    openvdb::Coord loc(x, y, z);
                    if (accessor.getValue(loc) > -offset)
                    {
                        yData.emplace_back(VOXELTOOL::VoxelEnumType::Outside);
                    }
                    else if (accessor.getValue(loc) == -background)
                    {
                        yData.emplace_back(VOXELTOOL::VoxelEnumType::Empty);
                    }
                    else
                    {
                        yData.emplace_back(VOXELTOOL::VoxelEnumType::Solid);
                    }
                }
                zData.emplace_back(yData);
            }
            voxelData.emplace_back(zData);
        }
    }

    void SelfSupportGenerator::topDown(const VOXELTOOL::VoxelBody &voxelBody, const int& x, const int& y, const int& z)
    {
        typename openvdb::FloatGrid::Accessor accessor = gridptr->getAccessor();
        // 向上进行寻找
        for (int zNow = z - 1; zNow < voxelData.size(); zNow--)
        {
            if (voxelBody.getVoxel(x, y, zNow) == VOXELTOOL::VoxelEnumType::SupportedSolid)
            {
                // 上下加自己体素往下掉
                for (int i = 1; i >= -1; i--)
                {
                    openvdb::Coord coordOrigin = toCoord(x, y, zNow + i);
                    ValueT valOrigin = accessor.getValue(coordOrigin);
                    openvdb::Coord coord = toCoord(x, y, z + i);
                    accessor.setValue(coord, valOrigin);
                }
                // 将自己之上的值设为背景值
                for (int z2 = z - 2; z2 >= zNow - 1; z2--)
                {
                    openvdb::Coord coord = toCoord(x, y, z2);
                    accessor.setValue(coord, background);
                }
                break;
            }
        }
    }

    void SelfSupportGenerator::setBackGround(const int& x, const int& y, const int& z, const VOXELTOOL::VoxelBody& voxelBody)
    {
        typename openvdb::FloatGrid::Accessor accessor = gridptr->getAccessor();
        // 往上搜索，如果搜索到支撑实体以后，将自己及以上到实体都设为背景值
        for (int zNow = z; zNow < voxelData.size(); zNow--)
        {
            if (voxelBody.getVoxel(x, y, zNow) == VOXELTOOL::VoxelEnumType::SupportedSolid)
            {
                // 将自己之上的值设为背景值
                for (int z2 = z; z2 >= zNow; z2--)
                {
                    setValue(x, y, z2, background);
                }
                break;
            }
        }
    }

    void SelfSupportGenerator::generate()
    {
        VOXELTOOL::VoxelBody voxelBody(voxelData, VOXELTOOL::PhysicsStrategy::Support5);
        voxelBody.support();
        typename openvdb::FloatGrid::Accessor accessor = gridptr->getAccessor();

        // 从上往下遍历
        int zRow = voxelData.size();
        int yRow = voxelData.at(0).size();
        int xRow = voxelData.at(0).at(0).size();
        for (int z = 0; z < zRow - 1; z++)
        {
            for (int y = 0; y < yRow; y++)
            {
                for (int x = 0; x < xRow; x++)
                {
                    if (voxelBody.getVoxel(x, y, z) == VOXELTOOL::VoxelEnumType::BorderSelfSupport)
                    {
                        // 如果下方是支撑实体，则全部变成背景值
                        if (voxelBody.getVoxel(x, y, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid ||
                            voxelBody.getVoxel(x - 1, y, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid ||
                            voxelBody.getVoxel(x, y - 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid ||
                            voxelBody.getVoxel(x + 1, y, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid ||
                            voxelBody.getVoxel(x, y + 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid/* ||
                            voxelBody.getVoxel(x - 1, y - 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid ||
                            voxelBody.getVoxel(x + 1, y - 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid ||
                            voxelBody.getVoxel(x + 1, y + 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid ||
                            voxelBody.getVoxel(x - 1, y + 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid*/
                            )
                        {
                            if (voxelBody.getVoxel(x, y, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid)
                            {
                                setBackGround(x, y, z, voxelBody);
                                setValue(x, y, z + 1, background);
                            }
                            else if (voxelBody.getVoxel(x - 1, y, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid)
                            {
                                setBackGround(x, y, z, voxelBody);
                                setValue(x - 1, y, z, background);
                                setValue(x - 1, y, z + 1, background);
                            }
                            else if (voxelBody.getVoxel(x, y - 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid)
                            {
                                setBackGround(x, y, z, voxelBody);
                                setValue(x, y - 1, z, background);
                                setValue(x, y - 1, z + 1, background);
                            }
                            else if (voxelBody.getVoxel(x + 1, y, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid)
                            {
                                setBackGround(x, y, z, voxelBody);
                                setValue(x + 1, y, z, background);
                                setValue(x + 1, y, z + 1, background);
                            }
                            else if (voxelBody.getVoxel(x, y + 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid)
                            {
                                setBackGround(x, y, z, voxelBody);
                                setValue(x, y + 1, z, background);
                                setValue(x, y + 1, z + 1, background);
                            }
                           /* else if (voxelBody.getVoxel(x - 1, y - 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid)
                            {
                                setBackGround(x, y, z, voxelBody);
                                setValue(x - 1, y - 1, z, background);
                                setValue(x - 1, y - 1, z + 1, background);
                            }
                            else if (voxelBody.getVoxel(x + 1, y - 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid)
                            {
                                setBackGround(x, y, z, voxelBody);
                                setValue(x + 1, y - 1, z, background);
                                setValue(x + 1, y - 1, z + 1, background);
                            }
                            else if (voxelBody.getVoxel(x + 1, y + 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid)
                            {
                                setBackGround(x, y, z, voxelBody);
                                setValue(x + 1, y + 1, z, background);
                                setValue(x + 1, y + 1, z + 1, background);
                            }
                            else if (voxelBody.getVoxel(x - 1, y + 1, z + 1) == VOXELTOOL::VoxelEnumType::SupportedSolid)
                            {
                                setBackGround(x, y, z, voxelBody);
                                setValue(x - 1, y + 1, z, background);
                                setValue(x - 1, y + 1, z + 1, background);
                            }*/
                        }
                        else
                        {
                            // 往上寻找
                            topDown(voxelBody, x, y, z);
                        }
                    }
                }
            }
        }

        // 把边界体素设成背景值
        for (int z = 0; z < zRow; z++)
        {
            for (int y = 0; y < yRow; y++)
            {
                for (int x = 0; x < xRow; x++)
                {
                    if (voxelBody.getVoxel(x, y, z) == VOXELTOOL::VoxelEnumType::SupportedSolid)
                    {
                        openvdb::Coord coord = toCoord(x, y, z);
                        ValueT valOrigin = accessor.getValue(coord);
                        accessor.setValue(coord, background);
                    }
                }
            }
        }
    }
}