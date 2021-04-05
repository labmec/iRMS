########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable(Paper2D main.cpp ${SOURCE} ${HEADERS})
target_include_directories(Paper2D PUBLIC ${common_dir})
# find_package(libInterpolate REQUIRED)
# find_package (Eigen3 3.3 REQUIRED)

find_package(BLAS REQUIRED)
target_link_libraries(Paper2D iMRS_library ${BLAS_LIBRARIES})
