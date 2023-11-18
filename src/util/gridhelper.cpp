#include "gridhelper.h"
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>

#include "tracer.h"

namespace ovdbutil
{    
    // TODO: Do I need to call initialize? Seems to work without it as well but the
    // docs say it should be called ones. It does a mutex lock-unlock sequence all
    // even if was called previously.
    FloatGridPtr mesh_to_grid(
        TracerInterrupter& tracer,
        trimesh::TriMesh* mesh,
        const openvdb::math::Transform& tr,
        float               exteriorBandWidth,
        float               interiorBandWidth,
        double voxel_size,
        int                 flags
    )
    {
        if (!mesh)
            return nullptr;

        openvdb::initialize();

        TriangleMeshDataAdapter adpater(*mesh);
        FloatGridPtr grid = openvdb::tools::meshToVolume<openvdb::FloatGrid>(
            tracer, adpater, tr, exteriorBandWidth,
            interiorBandWidth, voxel_size, flags);
    
        return grid;
    }
    
    FloatGridPtr mesh_to_grid(const Contour3D& mesh,
        const openvdb::math::Transform& tr,
        float exteriorBandWidth,
        float interiorBandWidth,
        int flags = 0)
    {
        openvdb::initialize();
        return openvdb::tools::meshToVolume<openvdb::FloatGrid>(
            Contour3DDataAdapter{ mesh }, tr, exteriorBandWidth, interiorBandWidth,
            flags);
    }
    
    Contour3D _volumeToMesh(FloatGridPtr grid,
        double      isovalue,
        double      adaptivity,
        bool        relaxDisorientedTriangles)
    {
        openvdb::initialize();
    
        Contour3D ret;

        if (grid)
        {
            std::vector<openvdb::Vec3s> points;
            std::vector<openvdb::Vec3I> triangles;
            std::vector<openvdb::Vec4I> quads;

            openvdb::tools::volumeToMesh(*grid, points, triangles, quads, isovalue,
                adaptivity, relaxDisorientedTriangles);


            ret.points.reserve(points.size());
            ret.faces3.reserve(triangles.size());
            ret.faces4.reserve(quads.size());

            for (auto& v : points) ret.points.emplace_back(to_vec3d(v));
            for (auto& v : triangles) ret.faces3.emplace_back(to_vec3i(v));
            for (auto& v : quads) ret.faces4.emplace_back(to_vec4i(v));
        }
    
        return ret;
    }
    
    trimesh::TriMesh* grid_to_mesh(FloatGridPtr grid,
        double                    isovalue,
        double                    adaptivity,
        bool                      relaxDisorientedTriangles)
    {
        return to_triangle_mesh(
            _volumeToMesh(grid, isovalue, adaptivity, relaxDisorientedTriangles));
    }
    
    Contour3D grid_to_contour3d(FloatGridPtr grid,
        double                    isovalue,
        double                    adaptivity,
        bool relaxDisorientedTriangles)
    {
        return _volumeToMesh(grid, isovalue, adaptivity,
            relaxDisorientedTriangles);
    }
}