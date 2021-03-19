########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable(MHM_Test main.cpp ${SOURCE} ${HEADERS} )
target_include_directories(MHM_Test PUBLIC ${common_dir})
# find_package(libInterpolate REQUIRED)
# find_package (Eigen3 3.3 REQUIRED)

target_link_libraries(MHM_Test iMRS_library blas)
