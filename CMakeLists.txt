set(OpenVDB_MAJOR_VERSION 9)
set(OpenVDB_MINOR_VERSION 0)
set(OpenVDB_PATCH_VERSION 1)
set(OpenVDB_VERSION "${OpenVDB_MAJOR_VERSION}.${OpenVDB_MINOR_VERSION}.${OpenVDB_PATCH_VERSION}")
set(OPENVDB_ABI_VERSION_NUMBER ${OpenVDB_MAJOR_VERSION})
set(CMAKE_CXX_STANDARD 14) 
include(GNUInstallDirs)

if(NOT OPENVDB_NAMESPACE_SUFFIX)
  set(OPENVDB_NAMESPACE_SUFFIX "" CACHE STRING
    "An optional custom string to append to the OpenVDB namespace. This is used to configure version.h.")
endif()

#########################################################################

message(STATUS "----------------------------------------------------")
message(STATUS "------------- Configuring OpenVDBCore --------------")
message(STATUS "----------------------------------------------------")

##########################################################################

# Collect and configure lib dependencies
__enable_bigobj()
#__assert_target(boost_iostreams)
__assert_target(tbb)

__cc_find(Trimesh2)
__assert_target(trimesh2)
__enable_bigobj()

# Set deps. Note that the order here is important. If we're building against
# Houdini 17.5 or later we must include IlmBase deps first to ensure the users
# chosen namespaced headers are correctly prioritized. Otherwise other include
# paths from shared installs (including houdini) may pull in the wrong headers

set(OPENVDB_CORE_DEPENDENT_LIBS "")

# We then choose to pull in TBB. If building aganst Houdini, TBB should
# always be pulled from there (we should always be using the version of TBB
# shipped with Houdini).

if(TARGET tbb)
    list(APPEND OPENVDB_CORE_DEPENDENT_LIBS tbb)
    add_definitions(-DUSE_TBB_MODULE)
endif()

set(USE_ZLIB 1)
if(USE_ZLIB)
  set(OPENVDB_USE_ZLIB 1)
  add_definitions(-DZ_PREFIX)
endif()

__assert_target(zlib)

list(APPEND OPENVDB_CORE_DEPENDENT_LIBS zlib)

if(NOT CC_BC_ANDROID)
	find_package(Threads REQUIRED)
	list(APPEND OPENVDB_CORE_DEPENDENT_LIBS Threads::Threads)
endif()

##### Core library configuration
# Configure VDB version in the form: ("%02x%02x%04x", major, minor, patch)
math(EXPR OPENVDB_PACKED_VERSION "${OpenVDB_MAJOR_VERSION} << 24")
math(EXPR OPENVDB_PACKED_VERSION "${OPENVDB_PACKED_VERSION} | ((${OpenVDB_MINOR_VERSION} & 0xFF) << 16)")
math(EXPR OPENVDB_PACKED_VERSION "${OPENVDB_PACKED_VERSION} | (${OpenVDB_PATCH_VERSION} & 0xFFFF)")

if(Imath_VERSION)
  set(OPENVDB_IMATH_VERSION ${Imath_VERSION})
endif()

configure_file(openvdb/version.h.in openvdb/version.h)

##########################################################################

set(OPENVDB_LIBRARY_SOURCE_FILES
  openvdb/Grid.cc
  openvdb/math/Half.cc
  openvdb/math/Maps.cc
  openvdb/math/Proximity.cc
  openvdb/math/Transform.cc
  openvdb/math/QuantizedUnitVec.cc
  openvdb/Metadata.cc
  openvdb/MetaMap.cc
  openvdb/openvdb.cc
  openvdb/io/Archive.cc
  openvdb/io/Compression.cc
  openvdb/io/GridDescriptor.cc
  openvdb/io/DelayedLoadMetadata.cc
  openvdb/Platform.cc
  openvdb/points/AttributeArray.cc
  openvdb/points/AttributeArrayString.cc
  openvdb/points/AttributeGroup.cc
  openvdb/points/AttributeSet.cc
  openvdb/points/points.cc
  openvdb/points/StreamCompression.cc
  openvdb/util/Formats.cc
  openvdb/util/Util.cc
)

if(NOT SHRINK_OPENVDB)
	list(APPEND OPENVDB_LIBRARY_SOURCE_FILES
		openvdb/Grid.cc
		openvdb/io/File.cc
		openvdb/io/Queue.cc
		openvdb/io/Stream.cc
		openvdb/io/TempFile.cc
		)
endif()

set(OPENVDB_LIBRARY_INCLUDE_FILES
  openvdb/Exceptions.h
  openvdb/Grid.h
  openvdb/Metadata.h
  openvdb/MetaMap.h
  openvdb/openvdb.h
  openvdb/Platform.h
  openvdb/PlatformConfig.h
  openvdb/Types.h
  openvdb/TypeList.h
)

set(OPENVDB_LIBRARY_IO_INCLUDE_FILES
  openvdb/io/Archive.h
  openvdb/io/Compression.h
  openvdb/io/DelayedLoadMetadata.h
  openvdb/io/File.h
  openvdb/io/GridDescriptor.h
  openvdb/io/io.h
  openvdb/io/Queue.h
  openvdb/io/Stream.h
  openvdb/io/TempFile.h
)

set(OPENVDB_LIBRARY_MATH_INCLUDE_FILES
  openvdb/math/BBox.h
  openvdb/math/ConjGradient.h
  openvdb/math/Coord.h
  openvdb/math/DDA.h
  openvdb/math/FiniteDifference.h
  openvdb/math/Half.h
  openvdb/math/HalfLimits.h
  openvdb/math/LegacyFrustum.h
  openvdb/math/Maps.h
  openvdb/math/Mat.h
  openvdb/math/Mat3.h
  openvdb/math/Mat4.h
  openvdb/math/Math.h
  openvdb/math/Operators.h
  openvdb/math/Proximity.h
  openvdb/math/QuantizedUnitVec.h
  openvdb/math/Quat.h
  openvdb/math/Ray.h
  openvdb/math/Stats.h
  openvdb/math/Stencils.h
  openvdb/math/Transform.h
  openvdb/math/Tuple.h
  openvdb/math/Vec2.h
  openvdb/math/Vec3.h
  openvdb/math/Vec4.h
)

set(OPENVDB_LIBRARY_POINTS_INCLUDE_FILES
  openvdb/points/AttributeArray.h
  openvdb/points/AttributeArrayString.h
  openvdb/points/AttributeGroup.h
  openvdb/points/AttributeSet.h
  openvdb/points/IndexFilter.h
  openvdb/points/IndexIterator.h
  openvdb/points/PointAdvect.h
  openvdb/points/PointAttribute.h
  openvdb/points/PointConversion.h
  openvdb/points/PointCount.h
  openvdb/points/PointDataGrid.h
  openvdb/points/PointDelete.h
  openvdb/points/PointGroup.h
  openvdb/points/PointMask.h
  openvdb/points/PointMove.h
  openvdb/points/PointRasterizeSDF.h
  openvdb/points/PointRasterizeTrilinear.h
  openvdb/points/PointSample.h
  openvdb/points/PointScatter.h
  openvdb/points/PointStatistics.h
  openvdb/points/PointTransfer.h
  openvdb/points/StreamCompression.h
)

set(OPENVDB_LIBRARY_TOOLS_INCLUDE_FILES
  openvdb/tools/Activate.h
  openvdb/tools/ChangeBackground.h
  openvdb/tools/Clip.h
  openvdb/tools/Composite.h
  openvdb/tools/Count.h
  openvdb/tools/Dense.h
  openvdb/tools/DenseSparseTools.h
  openvdb/tools/Diagnostics.h
  openvdb/tools/FastSweeping.h
  openvdb/tools/Filter.h
  openvdb/tools/FindActiveValues.h
  openvdb/tools/GridOperators.h
  openvdb/tools/GridTransformer.h
  openvdb/tools/Interpolation.h
  openvdb/tools/LevelSetAdvect.h
  openvdb/tools/LevelSetFilter.h
  openvdb/tools/LevelSetFracture.h
  openvdb/tools/LevelSetMeasure.h
  openvdb/tools/LevelSetMorph.h
  openvdb/tools/LevelSetPlatonic.h
  openvdb/tools/LevelSetRebuild.h
  openvdb/tools/LevelSetSphere.h
  openvdb/tools/LevelSetTracker.h
  openvdb/tools/LevelSetUtil.h
  openvdb/tools/Mask.h
  openvdb/tools/Merge.h
  openvdb/tools/MeshToVolume.h
  openvdb/tools/Morphology.h
  openvdb/tools/MultiResGrid.h
  openvdb/tools/NodeVisitor.h
  openvdb/tools/ParticleAtlas.h
  openvdb/tools/ParticlesToLevelSet.h
  openvdb/tools/PointAdvect.h
  openvdb/tools/PointIndexGrid.h
  openvdb/tools/PointPartitioner.h
  openvdb/tools/PointScatter.h
  openvdb/tools/PointsToMask.h
  openvdb/tools/PoissonSolver.h
  openvdb/tools/PotentialFlow.h
  openvdb/tools/Prune.h
  openvdb/tools/RayIntersector.h
  openvdb/tools/RayTracer.h
  openvdb/tools/SignedFloodFill.h
  openvdb/tools/Statistics.h
  openvdb/tools/TopologyToLevelSet.h
  openvdb/tools/ValueTransformer.h
  openvdb/tools/VectorTransformer.h
  openvdb/tools/VelocityFields.h
  openvdb/tools/VolumeAdvect.h
  openvdb/tools/VolumeToMesh.h
  openvdb/tools/VolumeToSpheres.h
)

set(OPENVDB_LIBRARY_TREE_INCLUDE_FILES
  openvdb/tree/InternalNode.h
  openvdb/tree/Iterator.h
  openvdb/tree/LeafBuffer.h
  openvdb/tree/LeafManager.h
  openvdb/tree/LeafNode.h
  openvdb/tree/LeafNodeBool.h
  openvdb/tree/LeafNodeMask.h
  openvdb/tree/NodeManager.h
  openvdb/tree/NodeUnion.h
  openvdb/tree/RootNode.h
  openvdb/tree/Tree.h
  openvdb/tree/TreeIterator.h
  openvdb/tree/ValueAccessor.h
)

set(OPENVDB_LIBRARY_UTIL_INCLUDE_FILES
  openvdb/util/CpuTimer.h
  openvdb/util/ExplicitInstantiation.h
  openvdb/util/Formats.h
  openvdb/util/logging.h
  openvdb/util/MapsUtil.h
  openvdb/util/Name.h
  openvdb/util/NodeMasks.h
  openvdb/util/NullInterrupter.h
  openvdb/util/PagedArray.h
  openvdb/util/Util.h
)

set(OPENVDB_LIBRARY_THREAD_INCLUDE_FILES
  openvdb/thread/Threading.h
)

##########################################################################

# Configure (static and shared lib) C flags

# Private defines (not inherited by dependent projects)
list(APPEND OPENVDB_CORE_PRIVATE_DEFINES OPENVDB_PRIVATE)

if(SHRINK_OPENVDB)
	list(APPEND OPENVDB_CORE_PRIVATE_DEFINES _STRINK_OPENVDB)
endif()

if(MINGW)
  list(APPEND OPENVDB_CORE_PRIVATE_DEFINES _USE_MATH_DEFINES)
endif()

set(OVDB_SRCS ${OPENVDB_LIBRARY_THREAD_INCLUDE_FILES}
		 ${OPENVDB_LIBRARY_UTIL_INCLUDE_FILES}
		 ${OPENVDB_LIBRARY_TREE_INCLUDE_FILES}
		 ${OPENVDB_LIBRARY_TOOLS_INCLUDE_FILES}
		 ${OPENVDB_LIBRARY_POINTS_INCLUDE_FILES}
		 ${OPENVDB_LIBRARY_MATH_INCLUDE_FILES}
		 ${OPENVDB_LIBRARY_IO_INCLUDE_FILES}
		 ${OPENVDB_LIBRARY_INCLUDE_FILES}
		 ${OPENVDB_LIBRARY_SOURCE_FILES}
		 )
set(INCS ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
set(DEFS ${OPENVDB_CORE_PRIVATE_DEFINES} OPENVDB_STATIC OPENVDB_STATICLIB)
set(LIBS ${OPENVDB_CORE_DEPENDENT_LIBS} boost_filesystem)

__add_real_target(openvdb lib SOURCE ${OVDB_SRCS}
							  INC ${INCS}
							  DEF ${DEFS}
							  LIB ${LIBS}
							  INTERFACE ${INCS}
							  )
			
			
set(SRCS include/ovdbutil/interface.h
		 include/ovdbutil/init.h
		 include/ovdbutil/hollowing.h
		 include/ovdbutil/subdivision.h
		 include/ovdbutil/booleanCom.h
		 src/init.cpp
		 src/hollowing.cpp
		 src/subdivision.cpp
		 src/booleanCom.cpp
		 src/util.h
		 src/util.cpp
		 )

if(HAVE_OVDBUTIL_BOOLEAN)
	list(APPEND SRCS include/ovdbutil/booleanCom.h
					 src/booleanCom.cpp
					)
endif()
		
list(APPEND INCS ${CMAKE_CURRENT_SOURCE_DIR}/include/)
list(APPEND LIBS trimesh2 mmesh boost_filesystem tbb openvdb)

__add_common_library(ovdbutil)

if(CC_INSTALL_OVDBUTIL)
	__install_directory_specif(${CMAKE_CURRENT_SOURCE_DIR}/include/ ovdbutil)
endif()

