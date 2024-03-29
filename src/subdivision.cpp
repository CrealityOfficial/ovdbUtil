#include "ovdbutil/subdivision.h"

#include "util.h"
#include "mmesh/trimesh/trimeshutil.h"
#include "ccglobal/tracer.h"
#include <openvdb/tools/RayIntersector.h>
#include "mmesh/create/createcylinder.h"
#include "mmesh/camera/ray.h"
#include "mmesh/util/adjacentoctree.h"
#include <openvdb/math/Vec3.h>
#include <openvdb/math/Coord.h>
#include "openvdb/tools/TopologyToLevelSet.h"

namespace ovdbutil
{
    trimesh::TriMesh* subdivison(trimesh::TriMesh* mesh,const SubdivisonParameter& parameter,const float sub, ccglobal::Tracer* tracer)
    {
        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(sub);
        std::vector<openvdb::math::Vec3s> s_points;
        std::vector<openvdb::math::Coord::Vec3I> s_faces;
        for (int i = 0; i < mesh->vertices.size(); i++)
        {
            s_points.push_back(openvdb::math::Vec3s(mesh->vertices.at(i).x /sub, mesh->vertices.at(i).y / sub, mesh->vertices.at(i).z / sub));
        }
        for (int i = 0; i < mesh->faces.size(); i++)
        {
            s_faces.push_back(openvdb::math::Coord::Vec3I(mesh->faces.at(i).x, mesh->faces.at(i).y, mesh->faces.at(i).z));
        }
        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_b(s_points, s_faces);
        //openvdb::FloatGrid::Ptr gridptr;
        openvdb::FloatGrid::Ptr subgrid = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_b, *xform);
         double iso_surface =0.;
         double adaptivity = 0.;
         auto omesh = grid_to_mesh(*subgrid, iso_surface, adaptivity, false);

        return omesh;
    }

    trimesh::TriMesh* remesh(trimesh::TriMesh* mesh, const float xformf, const int w, ccglobal::Tracer* tracer)
    {
        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(xformf);
        std::vector<openvdb::math::Vec3s> s_points;
        std::vector<openvdb::math::Coord::Vec3I> s_faces;
        for (int i = 0; i < mesh->vertices.size(); i++)
        {
            s_points.push_back(openvdb::math::Vec3s(mesh->vertices.at(i).x / xformf, mesh->vertices.at(i).y / xformf, mesh->vertices.at(i).z / xformf));
        }
        for (int i = 0; i < mesh->faces.size(); i++)
        {
            s_faces.push_back(openvdb::math::Coord::Vec3I(mesh->faces.at(i).x, mesh->faces.at(i).y, mesh->faces.at(i).z));
        }
        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_b(s_points, s_faces);
        //openvdb::FloatGrid::Ptr gridptr;
        openvdb::FloatGrid::Ptr subgrid = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_b, *xform);
        //openvdb::io::File("mypoints.vdb").write({ subgrid });
        subgrid = openvdb::tools::topologyToLevelSet(*subgrid, 3, 1, w, 3);   //extent    smooth steps
        double iso_surface = 0.;
        double adaptivity = 0.;
        auto omesh = grid_to_mesh(*subgrid, iso_surface, adaptivity, false);

        return omesh;
    }


}
