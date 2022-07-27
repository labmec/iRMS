#this allow us to define custom FindXXX.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# This is not the best. Please create (the extense) list of needed files
file(GLOB common_headers *.h json.hpp)
file(GLOB common_sources *.cpp)

# This package is here because iMRS is the only one that uses it
find_package(libInterpolate REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../libinterpolate_install ${CMAKE_SOURCE_DIR}/libinterpolate_install)

# install eigen
include(cmake/EnableEigen.cmake)
enable_eigen()

# Includes superLU using premaed .cmake files in folder Common/cmake
if(USING_SUPERLU)
    find_package(SuperLU REQUIRED)
endif(USING_SUPERLU)

# Create the iMRS library statically!
add_library(iMRS_library STATIC ${common_headers} ${common_sources})

#enabling TBB library
option(USING_TBB "Whether the TBB library will be linked in" OFF)
if(USING_TBB)
	find_package(Threads)
	find_package(TBB)
	target_compile_definitions(iMRS_library PRIVATE USING_TBB)
	target_link_libraries(iMRS_library PUBLIC TBB::tbb)
endif()

# Links the other libraries
target_link_libraries(iMRS_library PUBLIC NeoPZ::pz libInterpolate::Interpolate)
if(USING_SUPERLU)    
    target_link_libraries(iMRS_library PUBLIC ${SUPERLU_LIBRARIES})
endif(USING_SUPERLU)

# Links with dfn
target_link_libraries(iMRS_library PUBLIC DFN::dfnlib)
target_include_directories(iMRS_library PUBLIC ${DFN_INCLUDE_DIRS})

# Includes the directories
target_include_directories(iMRS_library PUBLIC ${PZ_INCLUDE_DIRS} 
    ${EIGEN3_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
if(USING_SUPERLU)    
    target_include_directories(iMRS_library PUBLIC ${SUPERLU_INCLUDES})
endif(USING_SUPERLU)


# create folder with meshes
set(RESERVOIRDIR ${CMAKE_SOURCE_DIR}/Reservoir)
set(FRACMESHESDIR ${CMAKE_SOURCE_DIR}/FracMeshes)
set (CONFIG_H "${CMAKE_SOURCE_DIR}/Common/imrs_config.h")
configure_file ("${CMAKE_SOURCE_DIR}/Common/imrs_config.h.in" ${CONFIG_H})
