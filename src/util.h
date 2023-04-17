#ifndef OVDBUTIL_UTIL_H
#define OVDBUTIL_UTIL_H
#include "trimesh2/TriMesh.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/GridTransformer.h>

namespace ccglobal
{
    class Tracer;
}

namespace ovdbutil
{
    /// @brief Contour3D
    struct Contour3D {
        std::vector<trimesh::vec3> points;
        std::vector<trimesh::ivec3> faces3;
        std::vector<trimesh::ivec4> faces4;

        Contour3D() = default;
        inline bool empty() const;
    };
    /// @brief Contour3DDataAdapter
    class Contour3DDataAdapter {
    public:
        const Contour3D& mesh;

        size_t polygonCount() const;
        size_t pointCount() const;
        size_t vertexCount(size_t n) const;

        // Return position pos in local grid index space for polygon n and vertex v
        void getIndexSpacePoint(size_t n, size_t v, openvdb::Vec3d& pos) const;
    };

    openvdb::FloatGrid::Ptr mesh_to_grid(trimesh::TriMesh* mesh,
        const openvdb::math::Transform& tr,
        float               exteriorBandWidth,
        float               interiorBandWidth,
        double voxel_size,
        int                 flags = 0,
        ccglobal::Tracer* tracer =nullptr
    );

    trimesh::TriMesh* grid_to_mesh(const openvdb::FloatGrid& grid,
        double                    isovalue,
        double                    adaptivity,
        bool                      relaxDisorientedTriangles);

    /// @brief TriangleMeshDataAdapter
    class TriangleMeshDataAdapter {
    public:
        trimesh::TriMesh& mesh;

        size_t polygonCount() const;
        size_t pointCount() const;
        size_t vertexCount(size_t) const;

        // Return position pos in local grid index space for polygon n and vertex v
        void getIndexSpacePoint(size_t n, size_t v, openvdb::Vec3d& pos) const;
    };




    trimesh::TriMesh* to_triangle_mesh(const Contour3D& ctour);

    template<class Grid>
    Contour3D _volumeToMesh(const Grid& grid,
        double      isovalue,
        double      adaptivity,
        bool        relaxDisorientedTriangles);

}

#endif // OVDBUTIL_UTIL_H