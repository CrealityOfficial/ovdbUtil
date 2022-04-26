__cc_find(Trimesh2)
__assert_target(trimesh2)
__enable_bigobj()

set(SRCS include/ovdbutil/interface.h
		 include/ovdbutil/init.h
		 include/ovdbutil/hollowing.h
		 src/init.cpp
		 src/hollowing.cpp
		 )
		
set(INCS ${CMAKE_CURRENT_SOURCE_DIR}/include/)
set(LIBS trimesh2 mmesh openvdb boost_iostreams tbb)

__add_common_library(ovdbutil)

if(CC_INSTALL_OVDBUTIL)
	__install_directory_specif(${CMAKE_CURRENT_SOURCE_DIR}/include/ ovdbutil)
endif()

