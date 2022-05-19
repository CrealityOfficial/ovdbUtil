#include "util.h"
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/GridTransformer.h>

#include "ccglobal/tracer.h"

namespace ovdbutil
{
    inline trimesh::vec3 to_vec3f(const openvdb::Vec3s& v) { return trimesh::vec3(v.x(), v.y(), v.z()); }
    inline trimesh::dvec3 to_vec3d(const openvdb::Vec3s& v) { return trimesh::dvec3(v.x(), v.y(), v.z()); }
    inline trimesh::ivec3 to_vec3i(const openvdb::Vec3I& v) { return trimesh::ivec3(int(v[0]), int(v[1]), int(v[2])); }
    inline trimesh::ivec4 to_vec4i(const openvdb::Vec4I& v) { return trimesh::ivec4(int(v[0]), int(v[1]), int(v[2]), int(v[3])); }
    
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
    
    // TODO: Do I need to call initialize? Seems to work without it as well but the
    // docs say it should be called ones. It does a mutex lock-unlock sequence all
    // even if was called previously.
    openvdb::FloatGrid::Ptr mesh_to_grid(trimesh::TriMesh* mesh,
        const openvdb::math::Transform& tr,
        float               exteriorBandWidth,
        float               interiorBandWidth,
        double voxel_size,
        int                 flags,
        ccglobal::Tracer* tracer
    )
    {
        openvdb::initialize();
    
        openvdb::FloatGrid::Ptr grid = openvdb::tools::meshToVolume<openvdb::FloatGrid>(
            TriangleMeshDataAdapter{ *mesh }, tr, exteriorBandWidth,
            interiorBandWidth, voxel_size, flags, tracer);
    
        return grid;
    }
    
    openvdb::FloatGrid::Ptr mesh_to_grid(const Contour3D& mesh,
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
    
    template<class Grid>
    Contour3D _volumeToMesh(const Grid& grid,
        double      isovalue,
        double      adaptivity,
        bool        relaxDisorientedTriangles)
    {
        openvdb::initialize();
    
        std::vector<openvdb::Vec3s> points;
        std::vector<openvdb::Vec3I> triangles;
        std::vector<openvdb::Vec4I> quads;
    
        openvdb::tools::volumeToMesh(grid, points, triangles, quads, isovalue,
            adaptivity, relaxDisorientedTriangles);
    
        Contour3D ret;
        ret.points.reserve(points.size());
        ret.faces3.reserve(triangles.size());
        ret.faces4.reserve(quads.size());
    
        for (auto& v : points) ret.points.emplace_back(to_vec3d(v));
        for (auto& v : triangles) ret.faces3.emplace_back(to_vec3i(v));
        for (auto& v : quads) ret.faces4.emplace_back(to_vec4i(v));
    
        return ret;
    }
    
    trimesh::TriMesh* grid_to_mesh(const openvdb::FloatGrid& grid,
        double                    isovalue,
        double                    adaptivity,
        bool                      relaxDisorientedTriangles)
    {
        return to_triangle_mesh(
            _volumeToMesh(grid, isovalue, adaptivity, relaxDisorientedTriangles));
    }
    
    Contour3D grid_to_contour3d(const openvdb::FloatGrid& grid,
        double                    isovalue,
        double                    adaptivity,
        bool relaxDisorientedTriangles)
    {
        return _volumeToMesh(grid, isovalue, adaptivity,
            relaxDisorientedTriangles);
    }
}