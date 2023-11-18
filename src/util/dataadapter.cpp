#include "dataadapter.h"

namespace ovdbutil
{    
    size_t Contour3DDataAdapter::polygonCount() const { return mesh.faces3.size() + mesh.faces4.size(); }
    size_t Contour3DDataAdapter::pointCount() const { return mesh.points.size(); }
    size_t Contour3DDataAdapter::vertexCount(size_t n) const { return n < mesh.faces3.size() ? 3 : 4; }
    
    
    inline bool Contour3D::empty() const
    {
        return points.empty() || (faces4.empty() && faces3.empty());
    }
    
    trimesh::TriMesh* to_triangle_mesh(const Contour3D& ctour) {
        trimesh::TriMesh* mesh = new trimesh::TriMesh();
        mesh->vertices = ctour.points;
        if (ctour.faces4.empty())
        {
            mesh->faces = ctour.faces3;
        }
        else
        {
            mesh->faces.reserve(ctour.faces3.size() + 2 * ctour.faces4.size());
            std::copy(ctour.faces3.begin(), ctour.faces3.end(),
                std::back_inserter(mesh->faces));
    
            for (const trimesh::ivec4& quad : ctour.faces4) {
                mesh->faces.emplace_back(quad[0], quad[1], quad[2]);
                mesh->faces.emplace_back(quad[2], quad[3], quad[0]);
            }
        }
    
        return mesh;
    }
    
    TriangleMeshDataAdapter::TriangleMeshDataAdapter(trimesh::TriMesh& m)
        :mesh(m)
    {

    }
    size_t TriangleMeshDataAdapter::polygonCount() const { return mesh.faces.size(); }
    size_t TriangleMeshDataAdapter::pointCount() const { return mesh.vertices.size(); }
    size_t TriangleMeshDataAdapter::vertexCount(size_t) const { return 3; }
    
    void TriangleMeshDataAdapter::getIndexSpacePoint(size_t          n,
        size_t          v,
        openvdb::Vec3d& pos) const
    {
        size_t vidx = size_t(mesh.faces[n][v]);
        trimesh::dvec3 p = trimesh::dvec3(mesh.vertices[vidx]);
        pos = { p.x, p.y, p.z };
    }
    
    void Contour3DDataAdapter::getIndexSpacePoint(size_t          n,
        size_t          v,
        openvdb::Vec3d& pos) const
    {
        size_t vidx = 0;
        if (n < mesh.faces3.size()) vidx = size_t(mesh.faces3[n][v]);
        else vidx = size_t(mesh.faces4[n - mesh.faces3.size()][v]);
    
        trimesh::vec3 p = mesh.points[vidx];
        pos = { p.x, p.y, p.z };
    }
}