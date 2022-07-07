#include "ovdbutil/booleanCom.h"
#include "openvdb/tools/TopologyToLevelSet.h"
#include "openvdb/tools/Morphology.h"
#include <openvdb/tools/Composite.h>
#include "util.h"
#include "mmesh/trimesh/trimeshutil.h"
#include "openvdb/tools/LevelSetUtil.h"
#include "ccglobal/tracer.h"
#include "mmesh/create/createcylinder.h"


namespace ovdbutil
{
    static trimesh::TriMesh* _voronia_boolcom(std::vector<trimesh::TriMesh>* mesh,
        const int type, const ovdbutil::BooleanParameter& param, ccglobal::Tracer* tracer)
    {
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.0f);

        //float D = 0;
        //float  out_range = param.externalWidth;
        //float  in_range = param.internalWidth;
        //float voxel_size = 0.05f;


        //openvdb::initialize();

        //// setup linear transform   
        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(0.25);

        openvdb::FloatGrid::Ptr grid;
        for (trimesh::TriMesh m : *mesh) {

            std::vector<openvdb::math::Vec3s> s_points;
            std::vector<openvdb::math::Coord::Vec3I> s_faces;
            for (int i = 0; i <m.vertices.size(); i++)
            {
                s_points.push_back(openvdb::math::Vec3s(m.vertices.at(i).x * 4, m.vertices.at(i).y * 4, m.vertices.at(i).z * 4));
            }
            for (int i = 0; i < m.faces.size(); i++)
            {
                s_faces.push_back(openvdb::math::Coord::Vec3I(m.faces.at(i).x, m.faces.at(i).y, m.faces.at(i).z));
            }


            openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_b(s_points, s_faces);
            //openvdb::FloatGrid::Ptr gridptr;
            openvdb::FloatGrid::Ptr subgrid = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_b, *xform);

            if (grid && subgrid) 
                openvdb::tools::csgUnion(*grid, *subgrid);
            else if (subgrid) 
                grid = std::move(subgrid);
        }
        //if (grid) {
        //    grid = openvdb::tools::levelSetRebuild(*grid, 0.);
        //}
        //else if (mesh.empty()) {
        //    // Splitting failed, fall back to hollow the original mesh
        //    grid = openvdb::tools::meshToVolume<openvdb::FloatGrid>(
        //        TriangleMeshDataAdapter{ mesh }, tr, exteriorBandWidth,
        //        interiorBandWidth, flags);
        //}

        if (tracer)
            tracer->progress(0.4f);
        if (!grid) {
            if (tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }

        if (tracer)
            tracer->progress(1.0f);
        double iso_surface = 0.;
        double adaptivity = 0.;
        auto omesh = ovdbutil::grid_to_mesh(*grid, iso_surface, adaptivity, false);

        return omesh;
    }


    static trimesh::TriMesh* _generate_boolcom(ovdbutil::TwoTrimesh* mesh,
        const int type, const ovdbutil::BooleanParameter& param, ccglobal::Tracer* tracer)
    {
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.0f);

        float D = 0;
        float  out_range = param.externalWidth;
        float  in_range = param.internalWidth;
        float voxel_size = 0.05f;
        

        openvdb::initialize();

        // setup linear transform   
        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(0.25);

        // CUBE 
        std::vector<openvdb::math::Vec3s> cube_points;
        std::vector<openvdb::math::Coord::Vec3I> cube_faces;
        for (int i = 0; i < mesh->m1->vertices.size(); i++)
        {
            cube_points.push_back(openvdb::math::Vec3s(mesh->m1->vertices.at(i).x * 4, mesh->m1->vertices.at(i).y * 4, mesh->m1->vertices.at(i).z * 4));
        }
        for (int i = 0; i < mesh->m1->faces.size(); i++)
        {
            cube_faces.push_back(openvdb::math::Coord::Vec3I(mesh->m1->faces.at(i).x, mesh->m1->faces.at(i).y, mesh->m1->faces.at(i).z));
        }

        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_a(cube_points, cube_faces);
        //openvdb::FloatGrid::Ptr gridptr;
        openvdb::FloatGrid::Ptr gridptr1 = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_a, *xform);

        // spherer 
        std::vector<openvdb::math::Vec3s> s_points;
        std::vector<openvdb::math::Coord::Vec3I> s_faces;
        for (int i = 0; i < mesh->m2->vertices.size(); i++)
        {
            s_points.push_back(openvdb::math::Vec3s(mesh->m2->vertices.at(i).x * 4, mesh->m2->vertices.at(i).y * 4, mesh->m2->vertices.at(i).z * 4));
        }
        for (int i = 0; i < mesh->m2->faces.size(); i++)
        {
            s_faces.push_back(openvdb::math::Coord::Vec3I(mesh->m2->faces.at(i).x, mesh->m2->faces.at(i).y, mesh->m2->faces.at(i).z));
        }

        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_b(s_points, s_faces);
        //openvdb::FloatGrid::Ptr gridptr;
        openvdb::FloatGrid::Ptr gridptr2 = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_b, *xform);
        //openvdb::FloatGrid::Ptr gridptrout1 = mesh_to_grid(mesh->m1, {}, out_range, in_range, voxel_size);



        if (type==0)  openvdb::tools::csgUnion(*gridptr1, *gridptr2);
        if (type == 1)  openvdb::tools::csgIntersection(*gridptr1, *gridptr2);
        if (type == 2)  openvdb::tools::csgDifference(*gridptr1, *gridptr2);
        if (type == 3)  openvdb::tools::csgDifference(*gridptr1, *gridptr2);

        if (tracer)
            tracer->progress(0.4f);
        if (!gridptr1) {
            if (tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }
        if (!gridptr2) {
            if (tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }

        if (tracer)
            tracer->progress(0.7f);
        


        if (tracer)
            tracer->progress(1.0f);
        double iso_surface = 0.;
        double adaptivity = 0.;
        auto omesh = ovdbutil:: grid_to_mesh(*gridptr1, iso_surface, adaptivity, false);

        return omesh;
    }
    trimesh::TriMesh* generateBoolcom( ovdbutil::TwoTrimesh* mesh,
        const int type,const ovdbutil::BooleanParameter& param, ccglobal::Tracer* tracer)
    {
        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;

//        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* meshptr = _generate_boolcom(mesh, type,param,tracer);

        if (meshptr) {
            mmesh::reverseTriMesh(meshptr);
        }

        return meshptr;
    }
    trimesh::TriMesh* voroniaBoolcomInputLines(std::vector<edge> &edges)
    {
        trimesh::TriMesh* outMesh = new trimesh::TriMesh();
        std::vector<trimesh::TriMesh> resuMerge;
        for (int i = 0; i < edges.size(); i++)
        {
            trimesh::vec3 centerPoint0(edges.at(i).p0);
            trimesh::vec3 centerPoint1(edges.at(i).p1);
            trimesh::vec3 vmiddle = (centerPoint0 + centerPoint1) / 2;
            trimesh::vec3 v0(centerPoint1 - centerPoint0);
            float len = trimesh::dist(centerPoint0, centerPoint1);
            len += 1.0;
            trimesh::TriMesh* symesh = mmesh:: createSoupCylinder(10, 1.6, len, vmiddle, v0);
            resuMerge.emplace_back(*symesh);
        }
        ovdbutil::BooleanParameter param;
        ccglobal::Tracer* tracer = nullptr;
        const int typeUnion = 0;
        return voroniaBoolcom(&resuMerge, 0, param, tracer);
    
    }
    trimesh::TriMesh* voroniaBoolcom(std::vector<trimesh::TriMesh>* mesh,
        const int type, const ovdbutil::BooleanParameter& param, ccglobal::Tracer* tracer)
    {
        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;

        //        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* meshptr = _voronia_boolcom(mesh, type, param, tracer);

        if (meshptr) {
            mmesh::reverseTriMesh(meshptr);
        }

        return meshptr;
    }

    static trimesh::TriMesh* _generate_Dilationcom(trimesh::TriMesh* mesh,
        const int& param, ccglobal::Tracer* tracer)
    {
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.0f);

        float offset = 1;
        float D = 0;
        float voxel_size = 0.05;


        openvdb::initialize();

        // setup linear transform   
        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform(0.25);

        // CUBE 
        std::vector<openvdb::math::Vec3s> cube_points;
        std::vector<openvdb::math::Coord::Vec3I> cube_faces;
        for (int i = 0; i < mesh->vertices.size(); i++)
        {
            cube_points.push_back(openvdb::math::Vec3s(mesh->vertices.at(i).x * 4, mesh->vertices.at(i).y * 4, mesh->vertices.at(i).z * 4));
        }
        for (int i = 0; i < mesh->faces.size(); i++)
        {
            cube_faces.push_back(openvdb::math::Coord::Vec3I(mesh->faces.at(i).x, mesh->faces.at(i).y, mesh->faces.at(i).z));
        }

        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_a(cube_points, cube_faces);
        //openvdb::FloatGrid::Ptr gridptr;
        openvdb::FloatGrid::Ptr gridptr1 = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_a, *xform);

        gridptr1 = openvdb::tools::topologyToLevelSet(*gridptr1, 3, 1, param, 1);   //change the  detail
        //openvdb::tools::dilateActiveValues(gridptr1.get()->tree(),param);

        if (tracer)
            tracer->progress(0.4f);
        if (!gridptr1) {
            if (tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }

        if (tracer)
            tracer->progress(0.7f);

        if (tracer)
            tracer->progress(1.0f);
        double iso_surface = 0.;
        double adaptivity = 0.;
        auto omesh = ovdbutil::grid_to_mesh(*gridptr1, iso_surface, adaptivity, false);

        return omesh;
    }

    trimesh::TriMesh* generateDilationcom(trimesh::TriMesh* mesh,
        const int param, ccglobal::Tracer* tracer)
    {
        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;

        //        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* meshptr =  _generate_Dilationcom(mesh,  param, tracer);

        //if (meshptr) {
        //    mmesh::reverseTriMesh(meshptr);
        //}

        return meshptr;
    }

    static trimesh::TriMesh* _generate_Erosioncom(trimesh::TriMesh* mesh,
        const int& param, ccglobal::Tracer* tracer)
    {
        if (tracer && tracer->interrupt())
            return nullptr;

        if (tracer)
            tracer->progress(0.0f);

        float offset = 1;
        float D = 0;
        float voxel_size = 0.05;


        openvdb::initialize();

        // setup linear transform   
        openvdb::math::Transform::Ptr xform = openvdb::math::Transform::createLinearTransform( 0.25);

        // CUBE 
        std::vector<openvdb::math::Vec3s> cube_points;
        std::vector<openvdb::math::Coord::Vec3I> cube_faces;
        for (int i = 0; i < mesh->vertices.size(); i++)
        {
            cube_points.push_back(openvdb::math::Vec3s(mesh->vertices.at(i).x * 4, mesh->vertices.at(i).y * 4, mesh->vertices.at(i).z *4));
        }
        for (int i = 0; i < mesh->faces.size(); i++)
        {
            cube_faces.push_back(openvdb::math::Coord::Vec3I(mesh->faces.at(i).x, mesh->faces.at(i).y, mesh->faces.at(i).z));
        }

        openvdb::tools::QuadAndTriangleDataAdapter<openvdb::math::Vec3s, openvdb::math::Coord::Vec3I> mesh_a(cube_points, cube_faces);
        //openvdb::FloatGrid::Ptr gridptr;
        openvdb::FloatGrid::Ptr gridptr1 = openvdb::tools::meshToVolume<openvdb::FloatGrid>(mesh_a, *xform);

      // openvdb:: tools::changeBackground(gridptr1->tree(), 1.0);

        if( 1) {
            auto maskPtr =openvdb:: tools::sdfInteriorMask(*gridptr1);
            gridptr1->topologyUnion(*maskPtr);
        }

       int  halfWidth = std::max(3, 1);
       int  closingSteps = std::max(1, 0);
       int  dilation = std::max(param, 0);

        openvdb:: MaskTree maskTree(gridptr1.get()->tree(), false/*background*/, openvdb::TopologyCopy());


        openvdb::tools::erodeActiveValues(maskTree, 3+param  , openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);
        openvdb::tools::dilateActiveValues(maskTree, closingSteps, openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);
        openvdb::tools::pruneInactive(maskTree);
        openvdb::tools::dilateActiveValues(maskTree, halfWidth, openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);
        openvdb::tools::erodeActiveValues(maskTree,  halfWidth, openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);
            
        openvdb::tools::pruneInactive(maskTree);


          //  M0.dilateVoxels(1+param,openvdb::tools:: NN_FACE, true);
          //  M0.erodeVoxels(1, openvdb::tools::NN_FACE, true);
        const float background = float(gridptr1->voxelSize()[0]) * float(halfWidth);
        typename openvdb::FloatTree::Ptr lsTree(
            new openvdb::FloatTree(maskTree, /*out=*/background, /*in=*/-background, openvdb::TopologyCopy()));


        //openvdb::tools::erodeActiveValues(maskTree, closingSteps, openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);
        //openvdb::tools::erodeActiveValues(maskTree, closingSteps, openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);
        //openvdb::tools::erodeActiveValues(*lsTree, closingSteps, openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);
        openvdb::tools::dilateActiveValues(*lsTree, 8+param, openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);
    // openvdb::tools::dilateActiveValues(*lsTree, closingSteps, openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);
        //openvdb::tools::dilateActiveValues(*lsTree, closingSteps, openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);

        lsTree->topologyDifference(maskTree);
        openvdb::tools::pruneLevelSet(*lsTree,  /*threading=*/true);

        // Create a level set grid from the tree
        typename openvdb::FloatGrid::Ptr lsGrid = openvdb::FloatGrid::create(lsTree);
        lsGrid->setTransform(gridptr1->transform().copy());
        lsGrid->setGridClass(openvdb::GRID_LEVEL_SET);

        // Use a PDE based scheme to propagate distance values from the
        // implicit zero crossing.
        openvdb::OPENVDB_VERSION_NAME::util::NullInterrupter* nullInterrupter = nullptr;
        openvdb::tools::ttls_internal::normalizeLevelSet(*lsGrid, 3 * halfWidth, nullInterrupter);

        // Additional filtering
        if (1 > 0) {
            openvdb::tools::ttls_internal::smoothLevelSet(*lsGrid, 1, halfWidth, nullInterrupter);
        }

        if (tracer)
            tracer->progress(0.4f);
        if (!gridptr1) {
            if (tracer)
                tracer->failed("Returned OpenVDB grid is NULL");
            return nullptr;
        }

        if (tracer)
            tracer->progress(0.7f);

        if (tracer)
            tracer->progress(1.0f);
        double iso_surface = 0;
        double adaptivity = 0.;
        auto omesh = ovdbutil::grid_to_mesh(*lsGrid, iso_surface, adaptivity, false);

        return omesh;
    }

    trimesh::TriMesh* generateErosioncom(trimesh::TriMesh* mesh,
        const int param, ccglobal::Tracer* tracer)
    {
        static const double MIN_OVERSAMPL = 3.;
        static const double MAX_OVERSAMPL = 8.;

        //        double voxel_scale = parameter.voxel_size_inout_range;
        trimesh::TriMesh* meshptr = _generate_Erosioncom(mesh, param, tracer);

        //if (meshptr) {
        //    mmesh::reverseTriMesh(meshptr);
        //}

        return meshptr;
    }
}