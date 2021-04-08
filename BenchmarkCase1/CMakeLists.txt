########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable(BenchmarkCase1 main.cpp ${SOURCE} ${HEADERS})
target_include_directories(BenchmarkCase1 PUBLIC ${common_dir})
# find_package(libInterpolate REQUIRED)
# find_package (Eigen3 3.3 REQUIRED)

find_package(BLAS)
target_link_libraries(BenchmarkCase1 iMRS_library ${BLAS_LIBRARIES})
