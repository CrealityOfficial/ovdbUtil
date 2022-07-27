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
namespace ovdbutil
{
    trimesh::TriMesh* subdivison(trimesh::TriMesh* mesh,const SubdivisonParameter& parameter, ccglobal::Tracer* tracer)
    {
        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(0.125);
        std::vector<openvdb::math::Vec3s> s_points;
        std::vector<openvdb::math::Coord::Vec3I> s_faces;
        for (int i = 0; i < mesh->vertices.size(); i++)
        {
            s_points.push_back(openvdb::math::Vec3s(mesh->vertices.at(i).x * 8, mesh->vertices.at(i).y * 8, mesh->vertices.at(i).z * 8));
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

}
