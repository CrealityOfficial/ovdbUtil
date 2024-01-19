#include "hollowing.h"
#include "ovdbutil/hollow.h"
#include "ovdbUtil/SelfSupportGenerator.h"

#include "util/gridhelper.h"
#include "util/tracer.h"

#include "mesh/trimeshutil.h"
#include "mesh/createcylinder.h"
#include "mesh/ray.h"
#include "mesh/adjacentoctree.h"
#include "msbase/mesh/dumplicate.h"
#include "msbase/mesh/get.h"
#include "../topomesh/internal/alg/volumeMesh.h"
#include "topomesh/interface/subdivision.h"
#include "topomesh/interface/letter.h"
#include "internal/alg/fillhoneycombs.h"
#include "topomesh/interface/utils.h"
#include "internal/data/mmesht.h"
#include "internal/data/CMesh.h"
#include "trimesh2/TriMesh_algo.h"
#include <openvdb/math/Vec3.h>
#include <openvdb/math/Coord.h>
#include <openvdb/tools/ValueTransformer.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/GridTransformer.h>
#include "VOXELTOOL/VoxelBody.hpp"
#include <limits>

#include <queue>

#include <iomanip>

namespace ovdbutil
{

    void FindTopZeroLevelSetAndGetBone(openvdb::FloatGrid::Ptr grid, std::vector<std::vector<bool>> &bone, int &c);

    trimesh::TriMesh *hollowMesh(trimesh::TriMesh *mesh,
                                 const HollowParameter &parameter, ccglobal::Tracer *tracer)
    {
        ovdbutil::HollowingParameter param;
        param.min_thickness = parameter.thickness;
        switch (parameter.precision)
        {
        case OutPrecesion::COARSE:
            param.precision = 1;
            break;
        case OutPrecesion::NORMAL:
            param.precision = 2;
            break;
        case OutPrecesion::ELABORATE:
            param.precision = 3;
            break;
        case OutPrecesion::EXTREME:
            param.precision = 4;
            break;
        }
        param.remain_main_shell = parameter.remain_main_shell;
        param.filter_tiny_shell = parameter.filter_tiny_shell;
        param.filter_shell = parameter.filter_shell;
        switch (parameter.style)
        {
        case HollowStyle::hs_none:
            param.fill_config.enable = false;
            break;
        case HollowStyle::hs_infill_grid:
            param.fill_config.enable = true;
            param.fill_config.fillratio = parameter.fill_ratio / 100.f;
            break;
            //-----empty----
        case HollowStyle::hs_self_support:
            param.fill_config.enable = false;
            param.self_support = true;
            break;
        }
        trimesh::TriMesh *hollowMesh = ovdbutil::hollowPrecisionMeshAndFill(mesh, param, tracer);
        return hollowMesh;
    }

    trimesh::TriMesh *shellMesh(trimesh::TriMesh *mesh, const ShellParameter &parameter, ccglobal::Tracer *tracer)
    {
        ovdbutil::HollowingParameter param;
        param.min_thickness = parameter.thickness;
        switch (parameter.precision)
        {
        case OutPrecesion::COARSE:
            param.precision = 1;
            break;
        case OutPrecesion::NORMAL:
            param.precision = 2;
            break;
        case OutPrecesion::ELABORATE:
            param.precision = 3;
            break;
        case OutPrecesion::EXTREME:
            param.precision = 4;
            break;
        }
        std::vector<int> faces = parameter.faces;
        trimesh::TriMesh *hollowMesh = SelectFacesHollow(mesh, faces, param, tracer);
        return hollowMesh;
    }

    openvdb::FloatGrid::Ptr redistance_grid(const openvdb::FloatGrid &grid, double iso, double er = 3.0, double ir = 3.0)
    {
        return openvdb::tools::levelSetRebuild(grid, float(iso), float(er), float(ir));
    }

    double Cross(const trimesh::vec2 &p1, const trimesh::vec2 &p2, const trimesh::vec2 &p3, const trimesh::vec2 &p4)
    {
        return (p2.x - p1.x) * (p4.y - p3.y) - (p2.y - p1.y) * (p4.x - p3.x);
    }

    double Area(const trimesh::vec2 &p1, const trimesh::vec2 &p2, const trimesh::vec2 &p3)
    {
        return Cross(p1, p2, p1, p3);
    }

    double fArea(const trimesh::vec2 &p1, const trimesh::vec2 &p2, const trimesh::vec2 &p3)
    {
        return fabs(Area(p1, p2, p3));
    }

    trimesh::vec2 Inter(const trimesh::vec2 &p1, const trimesh::vec2 &p2, const trimesh::vec2 &p3, const trimesh::vec2 &p4)
    {
        double k = fArea(p1, p2, p3) / fArea(p1, p2, p4);
        return trimesh::vec2((p3.x + k * p4.x) / (1 + k), (p3.y + k * p4.y) / (1 + k));
    }

    std::vector<mmesh::Ray> *generatorRay(trimesh::TriMesh *amesh, INNER_FILL_CONFIG fillConfig, const trimesh::vec3 &normal)
    {

        std::vector<mmesh::Ray> *result = new std::vector<mmesh::Ray>;
        std::vector<trimesh::vec2> Xlines;
        std::vector<trimesh::vec2> Ylines;

        amesh->need_bbox();
        float gap = fillConfig.fillRadius / fillConfig.fillratio;
        float minZ = amesh->bbox.min.z;
        float minX = amesh->bbox.min.x;
        float maxX = amesh->bbox.max.x;
        int startX = minX / gap + 1;
        int endX = maxX / gap;
        float minY = amesh->bbox.min.y;
        float maxY = amesh->bbox.max.y;
        int startY = minY / gap + 1;
        int endY = maxY / gap;

        // �ų������߽�Ľ���
        if (startX * gap - minX < fillConfig.fillRadius * 2)
        {
            startX++;
        }
        if (maxX - endX * gap < fillConfig.fillRadius * 2)
        {
            endX--;
        }
        if (startY * gap - minY < fillConfig.fillRadius * 2)
        {
            startY++;
        }
        if (maxY - endY * gap < fillConfig.fillRadius * 2)
        {
            endY--;
        }

        for (int n = startX; n <= endX; n++)
        {
            trimesh::vec2 pointStart(n * gap, minY);
            trimesh::vec2 pointEnd(n * gap, maxY);
            Xlines.push_back(pointStart);
            Xlines.push_back(pointEnd);
        }
        for (int n = startY; n <= endY; n++)
        {
            trimesh::vec2 pointStart(minX, n * gap);
            trimesh::vec2 pointEnd(maxX, n * gap);
            Ylines.push_back(pointStart);
            Ylines.push_back(pointEnd);
        }

        std::vector<trimesh::vec3> vctIntersection;
        for (int n = 1; n < Xlines.size(); n += 2)
        {
            // n-1,n
            for (int m = 1; m < Ylines.size(); m += 2)
            {
                trimesh::vec2 inters = Inter(Xlines[n - 1], Xlines[n], Ylines[m - 1], Ylines[m]);
                trimesh::point pointTemp(inters.x, inters.y, minZ);
                vctIntersection.push_back(pointTemp);
            }
        }

        for (trimesh::vec3 &apoint : vctIntersection)
        {
            mmesh::Ray aRay(apoint, normal);
            result->push_back(aRay);
        }
        return result;
    };

    static trimesh::TriMesh *_generate_interior(trimesh::TriMesh *mesh,
                                                double min_thickness,
                                                double voxel_scale,
                                                double closing_dist,
                                                ccglobal::Tracer *tracer, double voxel_size)
    {
        //_scale(voxel_scale, imesh);

        double offset = voxel_scale * min_thickness;
        double D = voxel_scale * closing_dist;
        float out_range = 0.03f * float(offset);
        float in_range = 1.9f * float(offset + D);

        TracerInterrupter interrupter(tracer);

        // bug for param flags: hollow Optimization
        openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(1.0);

        // std::cout << "trans size :" << transform->voxelSize() << "\n";
        // transform.postScale(0.5);
        // std::cout << "after trans size :" << transform->voxelSize() << "\n";

        std::vector<openvdb::math::Vec3s> cube_points;
        std::vector<openvdb::math::Coord::Vec3I> cube_faces;
        for (int i = 0; i < mesh->vertices.size(); i++)
        {
            cube_points.push_back(openvdb::math::Vec3s(mesh->vertices.at(i).x * 3, mesh->vertices.at(i).y * 3, mesh->vertices.at(i).z * 3));
        }
        for (int i = 0; i < mesh->faces.size(); i++)
        {
            cube_faces.push_back(openvdb::math::Coord::Vec3I(mesh->faces.at(i).x, mesh->faces.at(i).y, mesh->faces.at(i).z));
        }
        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_a(cube_points, cube_faces);

        trimesh::TriMesh *newmesh = new trimesh::TriMesh();
        *newmesh = *mesh;
        std::vector<int> selectfaces;
        // topomesh::findNeignborFacesOfSameAsNormal(newmesh, 80319, 1.f, selectfaces);

        std::vector<bool> delfaces(newmesh->faces.size(), false);
        std::vector<bool> markvexter(newmesh->vertices.size(), false);
        float zi;
        for (int fi : selectfaces)
        {
            delfaces[fi] = true;
            for (int vi = 0; vi < 3; vi++)
            {
                int v = newmesh->faces[fi][vi];
                if (!markvexter[v])
                {
                    markvexter[v] = true;
                    zi = newmesh->vertices[v].z;
                    newmesh->vertices[v].z = newmesh->vertices[v].z - 1.0f;
                }
            }
        }
        /* trimesh::remove_faces(mesh, delfaces);
         trimesh::remove_unused_vertices(mesh);*/

        // newmesh->write("vertex.ply");
        openvdb::FloatGrid::Ptr gridptr1 = mesh_to_grid(interrupter,
                                                        newmesh, *transform, 3.0, 3.0f, voxel_size, 0xE);
        // typename openvdb::FloatGrid::Accessor accessor1 = gridptr1->getAccessor();
        typename openvdb::FloatGrid::Accessor accessor1 = gridptr1->getAccessor();

        std::cout << "back ground : " << gridptr1->background() << "\n";

        /* openvdb::BoolGrid::Ptr boolptr = openvdb::BoolGrid::create();

         for (int fi = 0 ; fi < 500; fi++)
         {
             openvdb::Vec3s max(-200,-200,-200);
             openvdb::Vec3s min(200, 200, 200);

             for (int vi = 0; vi > 3; vi++)
             {
                 max.x() = std::max(max.x(), mesh->vertices[mesh->faces[fi][vi]].x);
                 max.y() = std::max(max.y(), mesh->vertices[mesh->faces[fi][vi]].y);
                 max.z() = std::max(max.z(), mesh->vertices[mesh->faces[fi][vi]].z);

                 min.x() = std::min(min.x(), mesh->vertices[mesh->faces[fi][vi]].x);
                 min.y() = std::min(min.y(), mesh->vertices[mesh->faces[fi][vi]].y);
                 min.z() = std::min(min.z(), mesh->vertices[mesh->faces[fi][vi]].z);
             }

             openvdb::CoordBBox bbx= openvdb::CoordBBox(min.x(),min.y(),min.z(), max.x(),max.y(),max.z());
             boolptr->fill(bbx,0.0);
         }
         openvdb::tools::VolumeToMesh mesher;*/
        // openvdb::tools::extractIsosurfaceMask();

        // openvdb::FloatGrid::Ptr gridptr1 = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_a, *transform, out_range, in_range, voxel_size, 0xE, tracer);
        /* openvdb::Mat4R xform;
         openvdb::tools::GridTransformer transformer(xform);
         transformer.transformGrid<openvdb::tools::PointSampler, openvdb::FloatGrid>(*gridptr1,);*/

        // openvdb::FloatGrid::Accessor accessor = gridptr1->getAccessor();
        // for (openvdb::FloatGrid::ValueOnCIter iter = gridptr1->cbeginValueOn(); iter; ++iter) {
        //     // Get the coordinates of the current voxel
        //     openvdb::Coord coord = iter.getCoord();
        //
        //     // Convert voxel coordinates to world coordinates
        //     openvdb::Vec3d worldPos = transform->indexToWorld(coord);
        //     std::cout << "Grid " << coord << " worldPos " << worldPos << std::endl;
        //     // Do something with the world position
        //     // ...
        // }
        // openvdb::CoordBBox bbx = gridptr1->evalActiveVoxelBoundingBox();
        // std::cout << " bbx.min().z() " << bbx.min().z() << " bbx.max().z() " << bbx.max().z() << std::endl;
        // gridptr1->clipGrid(openvdb::BBoxd(openvdb::math::Vec3i(bbx.min().x(),bbx.min().y(), -12), openvdb::math::Vec3i(bbx.max().x(), bbx.max().y(), bbx.max().z())));

        // openvdb::FloatGrid::Ptr gridptr =
        //     openvdb::tools::createLevelSetSphere<openvdb::FloatGrid>(
        //         /*radius=*/10.0, /*center=*/openvdb::Vec3f(0, 0, 0),
        //         /*voxel size=*/0.5, /*width=*/2.0);

        // openvdb::tools::signedFloodFillWithValues(gridptr->tree(), 2.0, -2.0);
        openvdb::Vec3f c = openvdb::Vec3f(0, 0, 0);
        openvdb::FloatGrid::Ptr grid =
            openvdb::FloatGrid::create(/*background value=*/2.0);
        using ValueT = typename openvdb::FloatGrid::ValueType;
        const ValueT outside = grid->background();
        const ValueT inside = -outside;
        int padding = int(openvdb::math::RoundUp(openvdb::math::Abs(outside)));
        std::cout << "background : " << grid->background() << "\n";
        std::cout << "padding : " << padding << "\n";

        std::cout << " grid->voxelSize() : " << grid->voxelSize() << "\n";
        int dim = int(10.f + padding);
        typename openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
        openvdb::Coord ijk;
        int &i = ijk[0], &j = ijk[1], &k = ijk[2];
        std::cout << "i : " << i << " j : " << j << " k : " << k << "\n";
        for (i = c[0] - dim; i < c[0] + dim; ++i)
        {
            const float x2 = openvdb::math::Pow2(i - c[0]);
            // const float x2 = openvdb::math::Pow(i - c[0],100);
            /*const float distx1 = openvdb::math::Pow2(i - 10.f);
            const float distx2 = openvdb::math::Pow2(i + 10.f);
            const float distx = std::min(distx1,distx2);*/
            for (j = c[1] - dim; j < c[1] + dim; ++j)
            {
                const float x2y2 = openvdb::math::Pow2(j - c[1]) + x2;
                // const float x2y2 = openvdb::math::Pow(j - c[1],100) + x2;
                /* const float disty1 = openvdb::math::Pow2(j - 10.f);
                 const float disty2 = openvdb::math::Pow2(j + 10.f);
                 const float disty = std::min(disty1, disty2);*/
                for (k = c[2] - dim; k < c[2] + dim; ++k)
                {
                    const float dist = openvdb::math::Sqrt(x2y2 + openvdb::math::Pow2(k - c[2])) - 10.f;
                    ValueT val = ValueT(dist);
                    if (val < inside || outside < val)
                        continue;
                    accessor.setValue(ijk, val);
                    /*  if (k == 0)
                      {
                          const float dist = openvdb::math::Sqrt(x2y2) - 3.f;
                          ValueT val = ValueT(dist);
                          if ( outside < val) continue;
                          accessor.setValue(ijk, val);
                      }*/
                    /* ValueT val = ValueT(i&j);
                     if (val < inside || outside < val) continue;
                     accessor.setValue(ijk, val);*/
                }
            }
        }

        openvdb::tools::signedFloodFill(grid->tree());
        /* for (openvdb::FloatGrid::ValueOnIter iter = grid->beginValueOn(); iter; ++iter) {

             openvdb::Coord coord = iter.getCoord();
             std::cout << "Grid " << coord << " value : " << iter.getValue() << std::endl;
         }*/
        ValueT backgound = grid->background();
        std::cout << "bbox : " << grid->evalActiveVoxelBoundingBox() << "\n";
        openvdb::Coord box1 = grid->evalActiveVoxelBoundingBox().getStart();
        openvdb::Coord box2 = grid->evalActiveVoxelBoundingBox().getEnd();

        std::vector<std::pair<int, int>> bone;
        int dis_x = box2.x() - box1.x() + 1;
        int dis_y = box2.y() - box1.y() + 1;
        int offset_x = std::abs(box1.x());
        int offset_y = std::abs(box1.y());
        std::vector<std::vector<bool>> mapp(dis_x, std::vector<bool>(dis_y, false));
        openvdb::Coord topc;
        int &x = topc[0], &y = topc[1], &z = topc[2];
        // z = (int)box2.z()-2;
        z = 0;
        //-----y----------
        for (x = box1.x(); x <= box2.x(); ++x)
        {
            int yy = 0, n = 0;
            for (y = box1.y(); y <= box2.y(); ++y)
            {
                // if (accessor.getValue(topc) != backgound)
                if (accessor.getValue(topc) <= 0.f)
                {
                    yy += y;
                    n++;
                    openvdb::Coord next(x, y + 1, z);
                    if (accessor.getValue(next) > 0.f)
                    {
                        int c = (int)(yy / n);
                        bone.push_back(std::make_pair(x, c));
                        mapp[x + offset_x][c + offset_y] = true;
                        yy = 0;
                        n = 0;
                    }
                }
            }
        }
        //-----------x------------
        for (y = box1.y(); y <= box2.y(); ++y)
        {
            int xx = 0, n = 0;
            for (x = box1.x(); x <= box2.y(); ++x)
            {
                // if (std::abs(accessor.getValue(topc)) != backgound)
                if (accessor.getValue(topc) <= 0.f)
                {
                    // std::cout << "* ";
                    std::cout << std::setprecision(3) << accessor.getValue(topc) << " ";
                    xx += x;
                    n++;
                    openvdb::Coord next(x + 1, y, z);
                    if (accessor.getValue(next) > 0.f)
                    {
                        int c = (int)(xx / n);
                        bone.push_back(std::make_pair(c, y));
                        mapp[c + offset_x][y + offset_y] = true;
                        xx = 0;
                        n = 0;
                    }
                }
                else
                {
                    // std::cout << "0 ";
                    std::cout << accessor.getValue(topc) << "    ";
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
        std::cout << "\n";
        for (int ii = box1.y(); ii <= box2.y(); ii++)
        {
            for (int jj = box1.x(); jj <= box2.x(); jj++)
            {
                bool prift = false;
                for (int v = 0; v < bone.size(); v++)
                {
                    if (bone[v].first == jj && bone[v].second == ii)
                    {
                        prift = true;
                        // std::cout << "* ";
                        break;
                    }
                }
                if (prift)
                    std::cout << "* ";
                else
                    std::cout << "0 ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
        std::cout << "\n";

        //-----use mark[][] save levelset
        int dz = z - 1;
        /* for (int dz = z - 1; dz >= box1.z(); dz--)
         {*/
        std::vector<std::pair<int, int>> mark;
        for (int dx = box1.x(); dx <= box2.x(); dx++)
        {
            for (int dy = box1.y(); dy <= box2.y(); dy++)
            {
                openvdb::Coord coord(dx, dy, dz);
                if (accessor.getValue(coord) < 0.f) // find zero
                {
                    // down step number
                    int step = 1;
                    bool is_break = false;
                    while (1)
                    {
                        // 回字 循环查找mark
                        for (int ddx = dx - step; ddx <= dx + step; ddx++)
                            for (int ddy = dy - step; ddy <= dy + step; ddy++)
                            {
                                if (std::sqrt(std::pow((std::abs(ddx - dx)), 2) + std::pow((std::abs(ddy - dy)), 2)) < step || ddx < box1.x() || ddx > box2.x() || ddy < box1.y() || ddy > box2.y())
                                    continue;
                                if (mapp[ddx + offset_x][ddy + offset_y])
                                {
                                    is_break = true;
                                    goto label;
                                }
                            }
                    label:
                        if (is_break)
                            break;
                        else
                        {
                            step++;
                        }
                        if (step >= dis_x / 2)
                        {
                            // find empty
                            break;
                        }
                    }
                    //---find z block  and  down
                }
                // if (!mapp[dx + offset_x][dy + offset_y]&& accessor.getValue(coord)!=backgound)
                //{
                //     bool prift = false;
                //     for (int ii = dx - 1; ii <= dx + 1; ii++)
                //     {
                //         for (int jj = dy - 1; jj <= dy + 1; jj++)
                //         {
                //             if (ii < box1.x() || jj < box1.y() || ii>box2.x() || jj>box2.y()) continue;
                //             if (mapp[ii + offset_x][jj + offset_y])
                //             {
                //                 //----isok------
                //                 prift = true;
                //                 break;
                //             }
                //         }
                //         if (prift)
                //             break;
                //     }
                //     if (prift)
                //     {
                //         std::cout << "* ";
                //         mark.push_back(std::make_pair(dx,dy));
                //     }
                //     else
                //     {
                //         std::cout << "0 ";
                //     }
                // }
                // else
                //{
                //     std::cout << "0 ";
                // }
            }
            std::cout << "\n";
        }
        for (int ii = 0; ii < mark.size(); ii++)
        {
            mapp[mark[ii].first + offset_x][mark[ii].second + offset_y] = true;
        }
        std::cout << "\n";
        std::cout << "\n";
        //  }

        /*openvdb::Coord ijk1;
        int& i1 = ijk1[0], & j1 = ijk1[1], & k1 = ijk1[2];
        for (k1 = box2.z(); k1 >box2.z()-1; --k1) {

            std::vector<std::pair<int, int>> map;
            for (i1 = box1.x(); i1 <= box2.x(); ++i1) {
                for (j1 = box1.y(); j1 <= box2.y(); ++j1) {
                    if (accessor.getValue(ijk1) != backgound)
                    {
                        int n = 0; int cn = 0;
                        for(int ii=i1-1;ii<=i1+1;ii++)
                            for (int jj = j1 - 1; jj <= j1 + 1; jj++)
                            {
                                if (ii >= box1.x() && ii <= box2.x() && jj >= box1.y() && jj <= box2.y())
                                {
                                    openvdb::Coord pad(ii,jj,k1);
                                    if (accessor.getValue(pad) != backgound)
                                    {
                                        n++;
                                        if (ii == i1 && (jj == j1 - 1 || jj == j1 + 1))
                                            cn++;
                                        if (jj == j1 && (ii == i1 - 1 || ii == i1 + 1))
                                            cn++;
                                    }
                                }
                            }
                        if (cn != 4&&n>=3)
                        {
                            map.push_back(std::make_pair(i1, j1));
                        }
                    }
                }
            }
            for (int p = 0; p < map.size(); p++)
            {
                if (k1 == box1.z()) continue;
                openvdb::Coord c(map[p].first,map[p].second,k1);
                openvdb::Coord b(map[p].first, map[p].second, k1-1);
                ValueT v = accessor.getValue(c);
                accessor.setValue(b,v);
                accessor.setValue(c, backgound);
            }
        }*/

        /* for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter) {
             float dist = iter.getValue();
             std::cout << "dist : " << dist <<" iter.getBoundingBox() :"<< iter.getBoundingBox()<<
                 " iter.getCoord() :"<< iter.getCoord()<<"iter.getLevel() : "<< iter.getLevel() <<"\n";

         }*/
        /* openvdb::CoordBBox cbbox(openvdb::Coord(0, 0, 0), openvdb::Coord(10, 10, 10));
         grid->sparseFill(cbbox, 0.);*/

        // openvdb::tools::meshToLevelSet();
        /* for (openvdb::FloatGrid::ValueOnIter iter = gridptr->beginValueOn(); iter; ++iter) {
             float dist = iter.getValue();
             iter.setValue((outside - dist) / width);
         }*/
        // gridptr->insertMeta("radius", openvdb::FloatMetadata(50.0));
        // std::cout << "Grid size :" << gridptr->voxelSize()<<" background : " <<gridptr->background();
        // gridptr->transformPtr()->postScale(0.5f);
        /* for (openvdb::FloatGrid::ValueOnIter iter = gridptr->beginValueOn(); iter; ++iter) {
             float dist = iter.getValue();
             iter.setValue(dist/2.0f);
         }*/
        // std::cout << "Grid size :" << gridptr->voxelSize() << " background : " << gridptr->background();
        /* for (openvdb::FloatGrid::ValueOnCIter iter = gridptr->cbeginValueOn(); iter; ++iter) {
             std::cout << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
         }*/
        // mesh->write("savemesh.ply");
        // openvdb::FloatGrid::Ptr gridptrout1 = mesh_to_grid(mesh, {}, out_range, in_range, voxel_size,0xF, tracer);
        //// Get the source and target grids' index space to world space transforms.
        // const openvdb::math::Transform
        //     & sourceXform = gridptr->transform(),
        //     & targetXform = gridptrout1->transform();
        //// Compute a source grid to target grid transform.
        //// (For this example, we assume that both grids' transforms are linear,
        //// so that they can be represented as 4 x 4 matrices.)
        // openvdb::Mat4R xform =
        //     sourceXform.createLinearTransform(40.0).get()->baseMap()->getAffineMap()->getMat4() *
        //     targetXform.createLinearTransform(40.0).get()->baseMap()->getAffineMap()->getMat4().inverse();
        // for (int i=0; i< 16; i++)
        //     xform.asPointer()[i]= xform.asPointer()[i]*10;
        //// Create the transformer.
        // openvdb::tools::GridTransformer transformer(xform);
        //// Resample using nearest-neighbor interpolation.
        // transformer.transformGrid<openvdb::tools::QuadraticSampler, openvdb::FloatGrid>(
        //     *gridptr, *gridptrout1);
        //// Prune the target tree for optimal sparsity.
        // gridptrout1->tree().prune();

        if (!gridptr1)
        {
            if (tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }

        if (tracer && tracer->interrupt())
            return nullptr;

        if (closing_dist > .0)
        {
            gridptr1 = redistance_grid(*gridptr1, -(offset + D), double(in_range));
        }
        else
        {
            D = -offset;
        }

        if (tracer)
        {
            tracer->progress(0.8f);
            if (tracer->interrupt())
            {
                return nullptr;
            }
        }

        // double iso_surface = D;
        double iso_surface = 0.;
        double adaptivity = 0.;
        auto omesh = grid_to_mesh(grid, iso_surface, adaptivity, true);
        omesh->write("omesh.ply");

        // openvdb::FloatGrid::Ptr gridptr2 = mesh_to_grid(omesh, *transform, out_range, in_range, voxel_size, 0xE, tracer);
        // openvdb::tools::csgDifference(*gridptr1, *gridptr2);
        ////_scale(1. / voxel_scale, omesh);
        // double iso_surface1 = 0.;
        // double adaptivity1 = 0.;
        // auto outmesh = ovdbutil::grid_to_mesh(*gridptr1, iso_surface1, adaptivity1, false);
        // outmesh->write("outmesh.ply");
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.95f);

        return omesh;
    }

    void FindTopZeroLevelSetAndGetBone(openvdb::FloatGrid::Ptr grid, std::vector<std::vector<bool>> &bone, int &c)
    {
        openvdb::Coord box_start = grid->evalActiveVoxelBoundingBox().getStart();
        openvdb::Coord box_end = grid->evalActiveVoxelBoundingBox().getEnd();

        typename openvdb::FloatGrid::Accessor accessor = grid->getAccessor();
        int offset_x = -box_start.x();
        int offset_y = -box_start.y();
        int dis_x = box_end.x() - box_start.x() + 1;
        int dis_y = box_end.y() - box_start.y() + 1;
        std::vector<std::vector<bool>> mark_container(dis_x, std::vector<bool>(dis_y, false));
        for (int z = box_end.z(); z >= box_start.z(); z--)
        {
            bool is_break = false;
            for (int y = box_start.y(); y <= box_end.y(); ++y)
            {
                for (int x = box_start.x(); x <= box_end.x(); ++x)
                {
                    if (accessor.getValue(openvdb::Coord(x, y, z)) <= 0.f && accessor.getValue(openvdb::Coord(x, y, z + 1)) > 0.f)
                    {
                        is_break = true;
                        c = z;
                        mark_container[x + offset_x][y + offset_y] = true;
                    }
                }
            }
            if (is_break)
                break;
        }

        // std::vector<std::vector<bool>> last_mark_container(dis_x, std::vector<bool>(dis_y, false));
        for (int i = 0; i < dis_y; i++)
        {
            std::vector<int> x_container;
            for (int j = 0; j < dis_x; j++)
            {
                if (mark_container[j][i])
                    x_container.push_back(j);
            }
            if (x_container.empty())
                continue;
            int x_center = 0;
            for (int xi = 0; xi < x_container.size(); xi++)
                x_center += x_container[xi];
            x_center /= x_container.size();
            bone[x_center][i] = true;
        }

        for (int i = 0; i < dis_x; i++)
        {
            std::vector<int> y_container;
            for (int j = 0; j < dis_y; j++)
            {
                if (mark_container[i][j])
                    y_container.push_back(j);
            }
            if (y_container.empty())
                continue;
            int y_center = 0;
            for (int yi = 0; yi < y_container.size(); yi++)
                y_center += y_container[yi];
            y_center /= y_container.size();
            bone[i][y_center] = true;
        }
    }

    trimesh::TriMesh *SelectFacesHollow(trimesh::TriMesh *mesh, std::vector<int> &selectfaces,
                                        const HollowingParameter &parameter, ccglobal::Tracer *tracer)
    {
        trimesh::point ave_normal(0, 0, 0);
        for (int fi : selectfaces)
        {
            trimesh::point normal = msbase::getFaceNormal(mesh, fi);
            ave_normal += normal;
        }
        ave_normal /= selectfaces.size();
        trimesh::normalize(ave_normal);

        trimesh::TriMesh *copymesh = new trimesh::TriMesh();
        *copymesh = *mesh;
        float extend = 2. * parameter.min_thickness;
        std::vector<bool> markvexter(copymesh->vertices.size(), false);
        std::vector<int> selectvexter;
        for (int fi : selectfaces)
        {
            for (int vi = 0; vi < 3; vi++)
            {
                int v = copymesh->faces[fi][vi];
                if (!markvexter[v])
                {
                    markvexter[v] = true;
                    copymesh->vertices[v] = copymesh->vertices[v] + extend * ave_normal;
                    selectvexter.push_back(v);
                }
            }
        }

        trimesh::fxform &&mat = trimesh::fxform::rot_into(ave_normal, trimesh::vec3(0, 0, -1));
        for (trimesh::point &p : mesh->vertices)
            p = mat * p;
        mesh->clear_bbox();
        mesh->need_bbox();
        float min_z = mesh->vertices[mesh->faces[selectfaces[0]][0]].z;

        mesh->need_neighbors();
        for (int vi = 0; vi < selectvexter.size(); vi++)
        {
            int v = selectvexter[vi];
            bool is_bound = false;
            for (int vv = 0; vv < mesh->neighbors[v].size(); vv++)
            {
                if (markvexter[mesh->neighbors[v][vv]])
                {
                    is_bound = true;
                    break;
                }
            }
            if (is_bound)
            {
                if (mesh->vertices[v].z < min_z)
                    min_z = mesh->vertices[v].z;
            }
        }

        for (trimesh::point &p : mesh->vertices)
            p.z = p.z - min_z;

        trimesh::TriMesh *returnmesh = hollowMesh(copymesh, parameter, tracer);
        copymesh->clear();
        if (!returnmesh || returnmesh->vertices.empty())
            return nullptr;
        std::vector<bool> delvertex(returnmesh->vertices.size(), false);
        for (trimesh::point &p : returnmesh->vertices)
        {
            p = mat * p;
            p.z = p.z - min_z;
        }

        trimesh::box2 box;
        box.min.x = std::numeric_limits<float>::max();
        box.min.y = std::numeric_limits<float>::max();
        box.max.x = std::numeric_limits<float>::min();
        box.max.y = std::numeric_limits<float>::min();
        for (int i : selectfaces)
        {
            for (int vi = 0; vi < 3; vi++)
            {
                float vx = mesh->vertices[mesh->faces[i][vi]].x;
                float vy = mesh->vertices[mesh->faces[i][vi]].y;
                if (vx < box.min.x)
                    box.min.x = vx;
                if (vx > box.max.x)
                    box.max.x = vx;
                if (vy < box.min.y)
                    box.min.y = vy;
                if (vy > box.max.y)
                    box.max.y = vy;
            }
        }

        for (int vi = 0; vi < returnmesh->vertices.size(); vi++)
        {
            trimesh::point v = returnmesh->vertices[vi];
            if (v.x < box.max.x && v.x > box.min.x && v.y < box.max.y && v.y > box.min.y && v.z > -1.2 * extend && v.z < 0.5f)
            {
                if (v.z < 0.f)
                    delvertex[vi] = true;
            }
        }
        trimesh::remove_vertices(returnmesh, delvertex);
        trimesh::remove_sliver_faces(returnmesh);

        std::vector<std::vector<trimesh::point>> sequentials = topomesh::GetOpenMeshBoundarys(returnmesh);
        std::vector<std::vector<trimesh::vec2>> totalpoly;
        for (int i = 0; i < sequentials.size(); i++)
        {
            totalpoly.push_back(std::vector<trimesh::vec2>());
            for (int j = 0; j < sequentials[i].size(); j++)
            {
                totalpoly[i].push_back(trimesh::vec2(sequentials[i][j].x, sequentials[i][j].y));
            }
        }
        trimesh::TriMesh *newmesh = new trimesh::TriMesh();
        *newmesh = *mesh;

        topomesh::embedingAndCutting(newmesh, totalpoly, selectfaces);
        std::vector<std::vector<std::vector<trimesh::vec2>>> Poly;
        Poly.push_back(totalpoly);
        std::vector<int> input;
        for (int i = 0; i < newmesh->faces.size(); i++)
            input.push_back(i);
        std::vector<int> output;
        topomesh::polygonInnerFaces(newmesh, Poly, input, output);
        std::vector<bool> innerfaces(newmesh->faces.size(), false);
        for (int i : output)
        {
            innerfaces[i] = true;
        }
        trimesh::remove_faces(newmesh, innerfaces);
        trimesh::remove_unused_vertices(newmesh);

        struct Equal_vec
        {
            bool operator()(const trimesh::vec2 &v1, const trimesh::vec2 &v2) const
            {
                return trimesh::len(v1 - v2) <= 1e-4;
            }
        };
        struct Hash_function
        {
            size_t operator()(const trimesh::vec2 &v) const
            {
                return (int(v.x * 1000) * 1000) + (int(v.y * 1000));
            }
        };
        int begin = selectvexter.size();
        std::unordered_map<trimesh::vec2, float, Hash_function, Equal_vec> refix;
        for (int vi = 0; vi < newmesh->vertices.size(); vi++)
        {
            trimesh::vec2 v = trimesh::vec2(newmesh->vertices[vi].x, newmesh->vertices[vi].y);
            refix.emplace(std::make_pair(v, newmesh->vertices[vi].z));
        }
        returnmesh->need_adjacentfaces();
        returnmesh->need_neighbors();

        for (int vi = 0; vi < returnmesh->vertices.size(); vi++)
        {
            if (returnmesh->neighbors[vi].size() != returnmesh->adjacentfaces[vi].size())
            {
                trimesh::vec2 v = trimesh::vec2(returnmesh->vertices[vi].x, returnmesh->vertices[vi].y);
                auto it = refix.find(v);
                if (it != refix.end())
                    returnmesh->vertices[vi].z = it->second;
            }
        }

        std::vector<bool> delfaces(mesh->faces.size(), false);
        for (int i : selectfaces)
            delfaces[i] = true;
        trimesh::remove_faces(mesh, delfaces);
        trimesh::remove_unused_vertices(mesh);

        int vertexsize = returnmesh->vertices.size();
        for (int vi = 0; vi < mesh->vertices.size(); vi++)
            returnmesh->vertices.push_back(mesh->vertices[vi]);
        for (int fi = 0; fi < mesh->faces.size(); fi++)
            returnmesh->faces.push_back(trimesh::TriMesh::Face(mesh->faces[fi][0] + vertexsize, mesh->faces[fi][1] + vertexsize, mesh->faces[fi][2] + vertexsize));

        vertexsize = returnmesh->vertices.size();
        for (int vi = 0; vi < newmesh->vertices.size(); vi++)
            returnmesh->vertices.push_back(newmesh->vertices[vi]);
        for (int fi = 0; fi < newmesh->faces.size(); fi++)
            returnmesh->faces.push_back(trimesh::TriMesh::Face(newmesh->faces[fi][0] + vertexsize, newmesh->faces[fi][1] + vertexsize, newmesh->faces[fi][2] + vertexsize));

        for (trimesh::point &p : returnmesh->vertices)
            p.z += min_z;
        trimesh::apply_xform(returnmesh, trimesh::xform::rot_into(trimesh::vec3(0, 0, -1), ave_normal));
        msbase::dumplicateMesh(returnmesh);
        msbase::mergeNearPoints(returnmesh, nullptr, 4e-3f);
        // returnmesh->write("resultmesh.ply");
        return returnmesh;
    }

    bool CheckConnectChunk(trimesh::TriMesh *mesh, std::vector<std::vector<int>> &chunks, std::vector<int> &block)
    {
        mesh->need_across_edge();
        std::vector<int> container;
        std::vector<bool> marked(mesh->faces.size(), false);
        for (int i = 0; i < chunks.size(); i++)
        {
            container.insert(container.end(), chunks[i].begin(), chunks[i].end());
            std::set<int> filter(container.begin(), container.end());
            container.assign(filter.begin(), filter.end());
        }
        for (int fi = 0; fi < container.size(); fi++)
            marked[container[fi]] = true;

        std::vector<int> temp;
        std::queue<int> que;
        que.push(container[0]);
        temp.push_back(container[0]);
        while (!que.empty())
        {
            marked[que.front()] = false;
            for (int i = 0; i < mesh->across_edge[que.front()].size(); i++)
            {
                int face = mesh->across_edge[que.front()][i];
                if (face == -1)
                    continue;
                if (marked[face])
                {
                    que.push(face);
                    marked[face] = false;
                    temp.push_back(face);
                }
            }
            que.pop();
        }
        if (temp.size() != container.size())
            return false;
        block.swap(container);
        return true;
    }

    trimesh::TriMesh *generateInterior(trimesh::TriMesh *mesh,
                                       const HollowingParameter &parameter, ccglobal::Tracer *tracer)
    {
        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;

        // I can't figure out how to increase the grid resolution through openvdb
        // API so the model will be scaled up before conversion and the result
        // scaled down. Voxels have a unit size. If I set voxelSize smaller, it
        // scales the whole geometry down, and doesn't increase the number of
        // voxels.
        //
        // max 8x upscale, min is native voxel size
        // double voxel_scale = MIN_OVERSAMPL + (MAX_OVERSAMPL - MIN_OVERSAMPL) * parameter.quality;
        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh *meshptr = _generate_interior(mesh, parameter.min_thickness, voxel_scale,
                                                       parameter.closing_distance, tracer, parameter.voxel_size);

        if (meshptr)
        {
            mmesh::reverseTriMesh(meshptr);
        }

        return meshptr;
    }

    OVDBUTIL_API trimesh::TriMesh *hollowMeshAndFill(trimesh::TriMesh *mesh,
                                                     const HollowingParameter &parameter, ccglobal::Tracer *tracer)
    {
        if (tracer)
        {
            tracer->progress(0.1f);
            if (tracer->interrupt())
            {
                return nullptr;
            }
        }

        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;
        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh *hollowMesh = _generate_interior(mesh, parameter.min_thickness, voxel_scale,
                                                          parameter.closing_distance, tracer, parameter.voxel_size);

        // ��Ǻ����0
        trimesh::TriMesh *outMesh = new trimesh::TriMesh();
        std::vector<trimesh::TriMesh *> meshtotalV;
        meshtotalV.push_back(mesh);
        if (hollowMesh != nullptr && hollowMesh->vertices.size() > 0)
        {
            meshtotalV.push_back(hollowMesh);
            if (1) //
            {
                trimesh::vec3 normal(0.5, 0.0, 1.0);
                std::vector<trimesh::TriMesh *> vctMesh = generateInfill(hollowMesh, normal, parameter);
                for (trimesh::TriMesh *cyMesh : vctMesh)
                {
                    meshtotalV.push_back(cyMesh);
                }
            }
            if (2) //
            {
                trimesh::vec3 normal(-0.5, 0.0, 1.0);
                std::vector<trimesh::TriMesh *> vctMesh = generateInfill(hollowMesh, normal, parameter);
                for (trimesh::TriMesh *cyMesh : vctMesh)
                {
                    meshtotalV.push_back(cyMesh);
                }
            }
            if (3) //
            {
                trimesh::vec3 normal(0.0, 0.5, 1.0);
                std::vector<trimesh::TriMesh *> vctMesh = generateInfill(hollowMesh, normal, parameter);
                for (trimesh::TriMesh *cyMesh : vctMesh)
                {
                    meshtotalV.push_back(cyMesh);
                }
            }
        }

        mmesh::mergeTriMesh(outMesh, meshtotalV);
        // outMesh->write("outmesh.ply");
        /*trimesh::TriMesh* newmesh = _generate_interior(outMesh, parameter.min_thickness, voxel_scale,
            parameter.closing_distance, tracer, parameter.voxel_size);
        newmesh->write("newmesh.ply");*/
        if (tracer)
        {
            tracer->progress(1.0f);
            if (tracer->interrupt())
            {
                return nullptr;
            }
        }
        return outMesh;
    }

    OVDBUTIL_API std::vector<trimesh::TriMesh *> generateInfill(trimesh::TriMesh *mesh, const trimesh::vec3 &normal, const HollowingParameter &param, ccglobal::Tracer *tracer /*= nullptr*/)
    {
        std::vector<trimesh::TriMesh *> vctMesh;
        if (!param.fill_config.enable)
        {
            return vctMesh;
        }

        std::vector<mmesh::Ray> *rays = generatorRay(mesh, param.fill_config, normal);
        ModelAdjacentOctree aModelAdj(mesh);
        for (mmesh::Ray &aray : *rays)
        {
            std::vector<std::pair<int, trimesh::dvec3>> faceIdIntersect;
            aModelAdj.lineCollide((trimesh::dvec3)aray.start, (trimesh::dvec3)trimesh::normalized(aray.dir), faceIdIntersect);

            // ������ͬ�ĵ�
            std::vector<trimesh::dvec3> vctIntersect;
            for (std::pair<int, trimesh::dvec3> &apair : faceIdIntersect)
            {
                static double preZ = apair.second.z - 1;
                if (preZ == apair.second.z)
                {
                    continue;
                }
                vctIntersect.push_back(apair.second);
                preZ = apair.second.z;
            }

            for (int n = 1; n < vctIntersect.size(); n += 2)
            {
                const trimesh::vec3 &startPoint = (trimesh::vec3)vctIntersect.at(n - 1);
                const trimesh::vec3 &endPoint = (trimesh::vec3)vctIntersect.at(n);
                float height = trimesh::distance(startPoint, endPoint) + param.min_thickness * 0.5; // param.min_thickness*0.25����������һ�㣬���õĸ���
                trimesh::vec3 centerPoint((endPoint.x + startPoint.x) * 0.5, (endPoint.y + startPoint.y) * 0.5, (endPoint.z + startPoint.z) * 0.5);
                trimesh::TriMesh *cylinderMesh = mmesh::createSoupCylinder(10, param.fill_config.fillRadius, height, centerPoint, aray.dir);
                vctMesh.push_back(cylinderMesh);
            }
        }
        return vctMesh;
    }

    trimesh::TriMesh *hollowPrecisionMeshAndFill(trimesh::TriMesh *mesh,
                                                 const HollowingParameter &parameter, ccglobal::Tracer *tracer)
    {
        if (tracer)
        {
            tracer->progress(0.1f);
            if (tracer->interrupt())
            {
                return nullptr;
            }
        }
        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;
        trimesh::TriMesh *returnmesh = hollowMesh(mesh, parameter, tracer);
        if (!returnmesh || returnmesh->vertices.empty())
            return nullptr;
        FindShellVolume(returnmesh, parameter);

        trimesh::TriMesh *outMesh = new trimesh::TriMesh();
        std::vector<trimesh::TriMesh *> meshtotalV;
        meshtotalV.push_back(mesh);
        if (returnmesh != nullptr && returnmesh->vertices.size() > 0)
        {
            meshtotalV.push_back(returnmesh);
            if (1) //
            {
                trimesh::vec3 normal(0.5, 0.0, 1.0);
                std::vector<trimesh::TriMesh *> vctMesh = generateInfill(returnmesh, normal, parameter);
                for (trimesh::TriMesh *cyMesh : vctMesh)
                {
                    meshtotalV.push_back(cyMesh);
                }
            }
            if (2) //
            {
                trimesh::vec3 normal(-0.5, 0.0, 1.0);
                std::vector<trimesh::TriMesh *> vctMesh = generateInfill(returnmesh, normal, parameter);
                for (trimesh::TriMesh *cyMesh : vctMesh)
                {
                    meshtotalV.push_back(cyMesh);
                }
            }
            if (3) //
            {
                trimesh::vec3 normal(0.0, 0.5, 1.0);
                std::vector<trimesh::TriMesh *> vctMesh = generateInfill(returnmesh, normal, parameter);
                for (trimesh::TriMesh *cyMesh : vctMesh)
                {
                    meshtotalV.push_back(cyMesh);
                }
            }
        }
        mmesh::mergeTriMesh(outMesh, meshtotalV);

        if (tracer)
        {
            tracer->progress(1.0f);
            if (tracer->interrupt())
            {
                return nullptr;
            }
        }
        return outMesh;
    }

    enum class SelfSupportStrategy : int
    {
        MidLine,       // 中轴线自支撑算法
        SelfAdaption   // 自适应自支撑算法
    };

    trimesh::TriMesh *hollowMesh(trimesh::TriMesh *mesh,
                                 const HollowingParameter &parameter, ccglobal::Tracer *tracer)
    {
        float min_thickness = parameter.min_thickness * parameter.precision * 1.0f;
        double voxel_scale = parameter.voxel_size_inout_range;
        double offset = voxel_scale * min_thickness;
        double D = voxel_scale * parameter.closing_distance;
        float out_range = 0.03f * float(offset);
        float in_range = 1.9f * float(offset + D);

        TracerInterrupter interrupter(tracer);

        openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(1.0f / (parameter.precision * 1.0f));
        std::vector<openvdb::math::Vec3s> cube_points;
        std::vector<openvdb::math::Coord::Vec3I> cube_faces;
        if (tracer)
            tracer->progress(0.2f);
        for (int i = 0; i < mesh->vertices.size(); i++)
        {
            cube_points.push_back(openvdb::math::Vec3s(mesh->vertices.at(i).x * parameter.precision, mesh->vertices.at(i).y * parameter.precision,
                                                       mesh->vertices.at(i).z * parameter.precision));
        }
        for (int i = 0; i < mesh->faces.size(); i++)
        {
            cube_faces.push_back(openvdb::math::Coord::Vec3I(mesh->faces.at(i).x, mesh->faces.at(i).y, mesh->faces.at(i).z));
        }
        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_a(cube_points, cube_faces);
        if (tracer)
            tracer->progress(0.4f);
        openvdb::FloatGrid::Ptr gridptr = openvdb::tools::meshToVolume<openvdb::FloatGrid>(interrupter,
                                                                                           mesh_a, *transform, out_range, in_range, parameter.voxel_size, 0xE);
        if (parameter.self_support)
        {
            SelfSupportStrategy strategy = SelfSupportStrategy::MidLine;
            switch (strategy)
            {
            case (SelfSupportStrategy::MidLine):
            {
                double inside_offset = -offset;
                using ValueT = typename openvdb::FloatGrid::ValueType;
                typename openvdb::FloatGrid::Accessor accessor = gridptr->getAccessor();
                ValueT backgound = gridptr->background();
                openvdb::Coord box1 = gridptr->evalActiveVoxelBoundingBox().getStart();
                openvdb::Coord box2 = gridptr->evalActiveVoxelBoundingBox().getEnd();
                int dis_x = box2.x() - box1.x() + 1;
                int dis_y = box2.y() - box1.y() + 1;
                int offset_x = -box1.x();
                int offset_y = -box1.y();
                int dis_z = box2.z() - box1.z() + 1;
                int offset_z = -box1.z();

                std::vector<std::vector<std::vector<bool>>> pendingcollose(dis_x, std::vector<std::vector<bool>>(dis_y, std::vector<bool>(dis_z, false)));
                std::vector<std::vector<std::vector<bool>>> mark_voxel(dis_x, std::vector<std::vector<bool>>(dis_y, std::vector<bool>(dis_z, false)));
                std::vector<std::vector<std::vector<bool>>> total_voxel(dis_x, std::vector<std::vector<bool>>(dis_y, std::vector<bool>(dis_z, false)));
                for (int zz = box2.z(); zz >= box1.z(); zz--)
                {
                    for (int yy = box1.y(); yy <= box2.y(); ++yy) {
                        for (int xx = box1.x(); xx <= box2.x(); ++xx) {
                            openvdb::Coord loc(xx, yy, zz);
                            if (accessor.getValue(loc) <= inside_offset)
                            {
                                pendingcollose[xx + offset_x][yy + offset_y][zz + offset_z] = true;
                                total_voxel[xx + offset_x][yy + offset_y][zz + offset_z] = true;
                            }
                            if (accessor.getValue(loc) <= inside_offset && (accessor.getValue(openvdb::Coord(xx, yy, zz + 1)) > inside_offset))
                            {
                                mark_voxel[xx + offset_x][yy + offset_y][zz + offset_z] = true;
                            }
                        }
                    }
                }

                std::vector<std::vector<std::vector<int>>> center_voxel(dis_x, std::vector<std::vector<int>>(dis_y, std::vector<int>(dis_z, 0)));
                for (int y = 0; y < pendingcollose[0].size(); y++)
                    for (int x = 0; x < pendingcollose.size(); x++)
                    {
                        std::vector<std::vector<int>> mark_z;
                        bool change = false;
                        for (int z = 0; z < pendingcollose[0][0].size(); z++)
                        {
                            if (!pendingcollose[x][y][z])
                            {
                                change = false;
                            }
                            else if (pendingcollose[x][y][z])
                            {
                                if (!change) mark_z.push_back(std::vector<int>());
                                change = true;
                                mark_z.back().push_back(z);
                            }
                        }

                        for (int mi = 0; mi < mark_z.size(); mi++)
                        {
                            int total_z = 0;
                            for (int mii = 0; mii < mark_z[mi].size(); mii++)
                                total_z += mark_z[mi][mii];
                            total_z /= mark_z[mi].size();
                            if (pendingcollose[x][y][total_z])
                                center_voxel[x][y][total_z]++;
                        }
                    }

                for (int z = 0; z < pendingcollose[0][0].size(); z++)
                    for (int y = 0; y < pendingcollose[0].size(); y++)
                    {
                        std::vector<std::vector<int>> mark_x;
                        bool change = false;
                        for (int x = 0; x < pendingcollose.size(); x++)
                        {
                            if (!pendingcollose[x][y][z])
                                change = false;
                            else if (pendingcollose[x][y][z])
                            {
                                if (!change) mark_x.push_back(std::vector<int>());
                                change = true;
                                mark_x.back().push_back(x);
                            }
                        }

                        for (int mi = 0; mi < mark_x.size(); mi++)
                        {
                            int total_x = 0;
                            for (int mii = 0; mii < mark_x[mi].size(); mii++)
                                total_x += mark_x[mi][mii];
                            total_x /= mark_x[mi].size();
                            if (pendingcollose[total_x][y][z])
                                center_voxel[total_x][y][z]++;
                        }
                    }

                for (int z = 0; z < pendingcollose[0][0].size(); z++)
                    for (int x = 0; x < pendingcollose.size(); x++)
                    {
                        std::vector<std::vector<int>> mark_y;
                        bool change = false;
                        for (int y = 0; y < pendingcollose[0].size(); y++)
                        {
                            if (!pendingcollose[x][y][z])
                                change = false;
                            else if (pendingcollose[x][y][z])
                            {
                                if (!change) mark_y.push_back(std::vector<int>());
                                change = true;
                                mark_y.back().push_back(y);
                            }
                        }

                        for (int mi = 0; mi < mark_y.size(); mi++)
                        {
                            int total_y = 0;
                            for (int mii = 0; mii < mark_y[mi].size(); mii++)
                                total_y += mark_y[mi][mii];
                            total_y /= mark_y[mi].size();
                            if (pendingcollose[x][total_y][z])
                                center_voxel[x][total_y][z]++;
                        }
                    }

                for (int z = 0; z < pendingcollose[0][0].size(); z++)
                    for (int y = 0; y < pendingcollose[0].size(); y++)
                        for (int x = 0; x < pendingcollose.size(); x++)
                            pendingcollose[x][y][z] = false;

                for (int z = 0; z < center_voxel[0][0].size(); z++)
                    for (int y = 0; y < center_voxel[0].size(); y++)
                        for (int x = 0; x < center_voxel.size(); x++)
                        {
                            if (center_voxel[x][y][z] >= 2)
                            {
                                int neighbor = 0;
                                for (int i = -1; i <= 1; i++)
                                    for (int j = -1; j <= 1; j++)
                                        for (int k = -1; k <= 1; k++)
                                        {
                                            if ((std::abs(i) + std::abs(j) + std::abs(k)) == 1)
                                            {
                                                if (x + i < 0 || y + j < 0 || z + k < 0) continue;
                                                if (center_voxel[x + i][y + j][z + k] >= 2)
                                                    neighbor++;
                                            }
                                        }
                                if (neighbor > 0)
                                {
                                    pendingcollose[x][y][z] = true;
                                }
                            }
                        }

                std::vector<std::vector<std::vector<bool>>> top_mark(dis_x, std::vector<std::vector<bool>>(dis_y, std::vector<bool>(dis_z, false)));
                for (int z = 0; z < pendingcollose[0][0].size(); z++)
                    for (int y = 0; y < pendingcollose[0].size(); y++)
                        for (int x = 0; x < pendingcollose.size(); x++)
                        {
                            if (pendingcollose[x][y][z])
                            {
                                for (int zz = z; zz < pendingcollose[0][0].size(); zz++)
                                {
                                    if (mark_voxel[x][y][zz] && !mark_voxel[x][y][zz - 1])
                                    {
                                        top_mark[x][y][zz] = true;
                                        break;
                                    }
                                }
                            }
                        }

                for (int z = top_mark[0][0].size() - 1; z >= 0; z--)
                    for (int y = 0; y < top_mark[0].size(); y++)
                        for (int x = 0; x < top_mark.size(); x++)
                        {
                            if (top_mark[x][y][z])
                            {
                                for (int i = -1; i <= 1; i++)
                                    for (int j = -1; j <= 1; j++)
                                    {
                                        if (total_voxel[x + i][y + j][z - 1])
                                            top_mark[x + i][y + j][z - 1] = true;
                                    }
                            }
                        }

                std::vector<std::vector<std::vector<bool>>> no_support(dis_x, std::vector<std::vector<bool>>(dis_y, std::vector<bool>(dis_z, false)));
                std::vector<std::vector<std::vector<bool>>> collose_shell(dis_x, std::vector<std::vector<bool>>(dis_y, std::vector<bool>(dis_z, false)));
                for (int z = top_mark[0][0].size() - 1; z >= 0; z--)
                    for (int y = 0; y < top_mark[0].size(); y++)
                        for (int x = 0; x < top_mark.size(); x++)
                        {
                            if (top_mark[x][y][z])
                            {
                                bool is_shell = false;
                                for (int i = -1; i <= 1; i++)
                                    for (int j = -1; j <= 1; j++)
                                        for (int k = -1; k <= 1; k++)
                                            if ((std::abs(i) + std::abs(j) + std::abs(k)) == 1)
                                                if (!top_mark[x + i][y + j][z + k])
                                                {
                                                    is_shell = true;
                                                }
                                if (is_shell)
                                    collose_shell[x][y][z] = true;
                            }
                        }

                for (int z = 0; z < collose_shell[0][0].size(); z++)
                    for (int y = 0; y < collose_shell[0].size(); y++)
                        for (int x = 0; x < collose_shell.size(); x++)
                        {
                            if (collose_shell[x][y][z])
                            {
                                bool is_support = false;
                                for (int i = -1; i <= 1; i++)
                                    for (int j = -1; j <= 1; j++)
                                    {
                                        if (collose_shell[x + i][y + j][z - 1])
                                            is_support = true;
                                        else
                                            if (!total_voxel[x + i][y + j][z - 1])
                                                is_support = true;
                                    }
                                if (!is_support)
                                    no_support[x][y][z] = true;
                            }
                        }

                std::vector<std::vector<std::vector<bool>>> no_support_2 = no_support;
                std::vector<std::vector<trimesh::ivec3>> line_support;
                for (int z = 0; z < no_support[0][0].size(); z++)
                    for (int y = 0; y < no_support[0].size(); y++)
                        for (int x = 0; x < no_support.size(); x++)
                        {
                            if (no_support[x][y][z])
                            {
                                line_support.push_back(std::vector<trimesh::ivec3>());
                                line_support.back().push_back(trimesh::ivec3(x, y, z));
                                no_support[x][y][z] = false;
                                std::queue<trimesh::ivec3> que;
                                que.push(trimesh::ivec3(x, y, z));
                                while (!que.empty())
                                {
                                    trimesh::ivec3 q = que.front();
                                    for (int i = -1; i <= 1; i++)
                                        for (int j = -1; j <= 1; j++)
                                        {
                                            if (q.x + i < 0 || q.y + j < 0) continue;
                                            if (i == 0 && j == 0) continue;
                                            if ((std::abs(i) + std::abs(j)) == 1)
                                                if (no_support[q.x + i][q.y + j][q.z])
                                                {
                                                    no_support[q.x + i][q.y + j][q.z] = false;
                                                    que.push(trimesh::ivec3(q.x + i, q.y + j, q.z));
                                                    line_support.back().push_back(trimesh::ivec3(q.x + i, q.y + j, q.z));
                                                }
                                        }
                                    que.pop();
                                }
                            }
                        }

                std::vector<trimesh::ivec3> volex_support;
                std::vector<std::vector<int>> lines_mark;
                for (int li = 0; li < line_support.size(); li++)
                {
                    bool is_intercross = false;
                    std::vector<int> line_mark;
                    for (int lii = 0; lii < line_support[li].size(); lii++)
                    {
                        trimesh::ivec3 p = line_support[li][lii];
                        int neighbor_mark = 0;
                        for (int i = -1; i <= 1; i++)
                            for (int j = -1; j <= 1; j++)
                            {
                                if (i == -1 && j == 0)
                                {
                                    if (no_support_2[p.x + i][p.y + j][p.z])
                                        neighbor_mark += 1;
                                }
                                else if (i == 0 && j == 1)
                                {
                                    if (no_support_2[p.x + i][p.y + j][p.z])
                                        neighbor_mark += 2;
                                }
                                else if (i == 1 && j == 0)
                                {
                                    if (no_support_2[p.x + i][p.y + j][p.z])
                                        neighbor_mark += 4;
                                }
                                else if (i == 0 && j == -1)
                                {
                                    if (no_support_2[p.x + i][p.y + j][p.z])
                                        neighbor_mark += 8;
                                }
                            }

                        if (neighbor_mark == 15)
                        {
                            is_intercross = true;
                            line_mark.push_back(2);

                            continue;
                        }
                        else if (neighbor_mark == 0)
                        {
                            line_mark.push_back(3);

                            continue;
                        }

                        int mf = neighbor_mark & 1;
                        int change_num = 0;
                        int mark_num = 0;
                        if (mf)
                        {
                            neighbor_mark += 16;
                        }
                        neighbor_mark >>= 1;
                        for (int i = 0; i < 4; i++)
                        {
                            int mk = neighbor_mark & 1;
                            if (mk)
                                mark_num++;
                            if (mk != mf)
                            {
                                mf = mk;
                                change_num++;
                            }
                            neighbor_mark >>= 1;
                        }


                        if (change_num == 2)
                        {
                            if (mark_num == 3)
                                line_mark.push_back(0);
                            else if (mark_num == 1 || mark_num == 2)
                                line_mark.push_back(1);
                            else
                                int aaa = 1;
                        }
                        else if (change_num == 4)
                            line_mark.push_back(0);

                    }
                    lines_mark.push_back(line_mark);
                    if (line_mark.size() == 1)
                    {
                        volex_support.push_back(trimesh::ivec3(line_support[li][0].x, line_support[li][0].y, line_support[li][0].z));
                        continue;
                    }

                    if (!is_intercross)
                    {
                        int cx = 0, cy = 0, dn = 0;
                        for (int ii = 0; ii < line_mark.size(); ii++)
                        {
                            if (line_mark[ii] == 1)
                            {
                                dn++;
                                cx += line_support[li][ii].x;
                                cy += line_support[li][ii].y;
                            }
                        }
                        cx /= dn;
                        cy /= dn;
                        volex_support.push_back(trimesh::ivec3(cx, cy, line_support[li][0].z));
                    }
                    else {
                        int cx = 0, cy = 0, dn = 0;
                        for (int ii = 0; ii < line_mark.size(); ii++)
                        {
                            if (line_mark[ii] == 2)
                            {
                                dn++;
                                cx += line_support[li][ii].x;
                                cy += line_support[li][ii].y;
                            }
                        }
                        cx /= dn;
                        cy /= dn;
                        volex_support.push_back(trimesh::ivec3(cx, cy, line_support[li][0].z));
                    }
                }

                for (int li = 0; li < line_support.size(); li++)
                {
                    for (int lii = 0; lii < line_support[li].size(); lii++)
                    {
                        if (lines_mark[li][lii] != 3 || lines_mark[li].size() != 1)
                        {
                            int dx = std::abs(line_support[li][lii].x - volex_support[li].x);
                            int dy = std::abs(line_support[li][lii].y - volex_support[li].y);
                            int dz = line_support[li][lii].z - dx - dy;
                            int tz = line_support[li][lii].z;
                            int zi = tz;
                            for (; zi >= dz; zi--)
                            {
                                if (!top_mark[line_support[li][lii].x][line_support[li][lii].y][zi])
                                    break;
                                top_mark[line_support[li][lii].x][line_support[li][lii].y][zi] = false;
                            }
                            if (zi == 0) continue;
                            if (lines_mark[li][lii] == 1 && top_mark[line_support[li][lii].x][line_support[li][lii].y][zi])
                            {
                                for (int zii = zi; zii >= 0; zii--)
                                {
                                    if (!top_mark[line_support[li][lii].x][line_support[li][lii].y][zii])
                                        break;
                                    top_mark[line_support[li][lii].x][line_support[li][lii].y][zii] = false;
                                }
                            }
                            if (lines_mark[li][lii] == 1)
                            {
                                std::vector<trimesh::ivec3> volex_neighbor;
                                volex_neighbor.push_back(trimesh::ivec3(line_support[li][lii].x - 1, line_support[li][lii].y, line_support[li][lii].z));
                                volex_neighbor.push_back(trimesh::ivec3(line_support[li][lii].x + 1, line_support[li][lii].y, line_support[li][lii].z));
                                volex_neighbor.push_back(trimesh::ivec3(line_support[li][lii].x, line_support[li][lii].y - 1, line_support[li][lii].z));
                                volex_neighbor.push_back(trimesh::ivec3(line_support[li][lii].x, line_support[li][lii].y + 1, line_support[li][lii].z));
                                for (int neighbor_i = 0; neighbor_i < 4; neighbor_i++)
                                {
                                    const trimesh::ivec3& p = volex_neighbor[neighbor_i];
                                    if (p.x < 0 || p.x >= top_mark.size() || p.y < 0 || p.y >= top_mark[0].size()) continue;
                                    int zu = p.z + 1;
                                    int zd = p.z - 1;
                                    top_mark[p.x][p.y][p.z] = false;
                                    while (1)
                                    {
                                        if (!top_mark[p.x][p.y][zu])
                                            break;
                                        top_mark[p.x][p.y][zu] = false;
                                        zu++;
                                    }
                                    while (1)
                                    {
                                        if (!top_mark[p.x][p.y][zd])
                                            break;
                                        top_mark[p.x][p.y][zd] = false;
                                        zd--;
                                    }
                                }
                            }
                        }
                    }
                }

                for (int z = box2.z(); z >= box1.z(); z--)
                    for (int y = box1.y(); y <= box2.y(); ++y) {
                        for (int x = box1.x(); x <= box2.x(); ++x) {
                            openvdb::Coord loc(x, y, z);
                            openvdb::Coord up(x, y, z + 1);
                            //openvdb::Coord down(x, y, z - 1);
                            ValueT val_loc = accessor.getValue(loc);
                            ValueT val_up = accessor.getValue(up);
                            // ValueT val_down = accessor.getValue(down);
                            if (accessor.getValue(loc) <= inside_offset && accessor.getValue(up) > inside_offset) {
                                if (!top_mark[x + offset_x][y + offset_y][z + offset_z])
                                {
                                    int zi = z + offset_z;
                                    for (; zi >= 0; zi--)
                                    {
                                        if (accessor.getValue(openvdb::Coord(x, y, zi - offset_z)) > inside_offset)
                                            break;
                                        if (top_mark[x + offset_x][y + offset_y][zi])
                                        {
                                            accessor.setValue(openvdb::Coord(x, y, zi - offset_z + 1), val_up);
                                            accessor.setValue(openvdb::Coord(x, y, zi - offset_z), val_loc);
                                            //accessor.setValue(openvdb::Coord(x, y, zi - offset_z-1), val_down);
                                            break;
                                        }
                                        else
                                            accessor.setValue(openvdb::Coord(x, y, zi - offset_z), backgound);
                                    }
                                }
                            }
                        }
                    }
                break;
            }
            case (SelfSupportStrategy::SelfAdaption):
            {
                // 自适应自支撑算法
                SelfSupportGenerator generator(gridptr, offset);
                generator.toVoxelEnumType();
                generator.generate();
                break;
            }
            }
        }

        if (!gridptr)
        {
            if (tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }

        if (tracer && tracer->interrupt())
            return nullptr;

        if (parameter.closing_distance > .0)
        {
            gridptr = redistance_grid(*gridptr, -(offset + D), double(in_range));
        }
        else
        {
            D = -offset;
        }

        double iso_surface = D;
        double adaptivity = 0.;
        if (tracer)
            tracer->progress(0.8f);

        trimesh::TriMesh *hollowMesh = grid_to_mesh(gridptr, iso_surface, adaptivity, false);
        hollowMesh->write("realhollowMesh.ply");
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.9f);

        return hollowMesh;
    }

    void FindShellVolume(trimesh::TriMesh *mesh, const HollowingParameter &parameter)
    {
        if (!mesh)
            return;
        mesh->clear_across_edge();
        mesh->need_across_edge();
        std::vector<bool> mark(mesh->faces.size(), false);
        std::vector<std::vector<int>> vol_container;
        int p = 0;
        bool pass = true;
        while (pass)
        {
            std::queue<int> facequeue;
            facequeue.push(p);
            std::vector<int> result;
            while (!facequeue.empty())
            {
                if (facequeue.front() == -1 || mark[facequeue.front()])
                {
                    facequeue.pop();
                    continue;
                }
                mark[facequeue.front()] = true;
                result.push_back(facequeue.front());
                for (int i = 0; i < mesh->across_edge[facequeue.front()].size(); i++)
                {
                    int face = mesh->across_edge[facequeue.front()][i];
                    if (face == -1 || mark[face])
                        continue;
                    facequeue.push(face);
                }
                facequeue.pop();
            }
            /* float vol = topomesh::getMeshVolume(mesh, result);
             vol = std::abs(vol);
             if (vol <= volume)
             {
                 del.insert(del.end(),result.begin(),result.end());
             }*/
            vol_container.push_back(result);
            int fi = 0;
            for (; fi < mesh->faces.size(); fi++)
            {
                if (!mark[fi])
                {
                    p = fi;
                    break;
                }
            }
            if (fi == mesh->faces.size())
                pass = false;
        }

        if (parameter.remain_main_shell)
        {
            std::sort(vol_container.begin(), vol_container.end(), [&](std::vector<int> &a, std::vector<int> &b) -> bool
                      {
                    float va = topomesh::getMeshVolume(mesh, a);
                    float vb = topomesh::getMeshVolume(mesh, b);
                    return std::abs(va) > std::abs(vb); });
            std::vector<bool> delfaces(mesh->faces.size(), false);
            for (int c = 1; c < vol_container.size(); c++)
            {
                for (int i = 0; i < vol_container[c].size(); i++)
                    delfaces[vol_container[c][i]] = true;
            }
            trimesh::remove_faces(mesh, delfaces);
            trimesh::remove_unused_vertices(mesh);
        }
        else if (parameter.filter_shell)
        {
            std::vector<bool> delfaces(mesh->faces.size(), false);
            for (int c = 0; c < vol_container.size(); c++)
            {
                float v = topomesh::getMeshVolume(mesh, vol_container[c]);
                v = std::abs(v);
                if (v <= parameter.filter_tiny_shell)
                    for (int i = 0; i < vol_container[c].size(); i++)
                        delfaces[vol_container[c][i]] = true;
            }
            trimesh::remove_faces(mesh, delfaces);
            trimesh::remove_unused_vertices(mesh);
        }
    }
}
