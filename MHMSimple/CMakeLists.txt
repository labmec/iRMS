########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable(MHMSimple main.cpp ${SOURCE} ${HEADERS} )
target_include_directories(MHMSimple PUBLIC ${common_dir})
# find_package(libInterpolate REQUIRED)
# find_package (Eigen3 3.3 REQUIRED)

target_link_libraries(MHMSimple iMRS_library blas)
