#ifndef DATAADAPTER_OVDBUTIL_UTIL_H
#define DATAADAPTER_OVDBUTIL_UTIL_H
#include "trimesh2/TriMesh.h"

#include <openvdb/Types.h>

namespace ovdbutil
{
    inline trimesh::vec3 to_vec3f(const openvdb::Vec3s& v) { return trimesh::vec3(v.x(), v.y(), v.z()); }
    inline trimesh::dvec3 to_vec3d(const openvdb::Vec3s& v) { return trimesh::dvec3(v.x(), v.y(), v.z()); }
    inline trimesh::ivec3 to_vec3i(const openvdb::Vec3I& v) { return trimesh::ivec3(int(v[0]), int(v[1]), int(v[2])); }
    inline trimesh::ivec4 to_vec4i(const openvdb::Vec4I& v) { return trimesh::ivec4(int(v[0]), int(v[1]), int(v[2]), int(v[3])); }

    /// @brief Contour3D
    struct Contour3D {
        std::vector<trimesh::vec3> points;
        std::vector<trimesh::ivec3> faces3;
        std::vector<trimesh::ivec4> faces4;

        Contour3D() = default;
        inline bool empty() const;
    };

    trimesh::TriMesh* to_triangle_mesh(const Contour3D& ctour);

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

    /// @brief TriangleMeshDataAdapter
    class TriangleMeshDataAdapter {
    public:
        TriangleMeshDataAdapter(trimesh::TriMesh& m);
        ~TriangleMeshDataAdapter() {}

        size_t polygonCount() const;
        size_t pointCount() const;
        size_t vertexCount(size_t) const;

        // Return position pos in local grid index space for polygon n and vertex v
        void getIndexSpacePoint(size_t n, size_t v, openvdb::Vec3d& pos) const;
    protected:
        trimesh::TriMesh& mesh;
    };
}

#endif // DATAADAPTER_OVDBUTIL_UTIL_H