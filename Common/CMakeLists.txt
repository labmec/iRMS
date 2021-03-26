########### next target ###############


file(GLOB common_headers *.h)
file(GLOB common_sources *.cpp)



set(common_files ${common_headers} ${common_sources} PARENT_SCOPE)
set(common_dir ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

find_package(libInterpolate REQUIRED)
find_package (Eigen3 3.3 REQUIRED)


add_library(iMRS_library ${common_headers} ${common_sources})
target_link_libraries(iMRS_library NeoPZ::pz libInterpolate::Interpolate)

if(USING_SUPERLU)
    include(SuperLU.cmake)
    target_link_libraries(iMRS_library ${SUPERLU_LIBRARIES})
endif(USING_SUPERLU)
# if(USING_SUPERLU)
#     include(SuperLU.cmake)
#     target_link_libraries(iMRS_library NeoPZ::pz libInterpolate::Interpolate ${SUPERLU_LIBRARIES})
# else(USING_SUPERLU)
#     target_link_libraries(iMRS_library NeoPZ::pz libInterpolate::Interpolate)
# endif(USING_SUPERLU)
