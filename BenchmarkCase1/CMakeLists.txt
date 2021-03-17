########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable(BecnhmarkCase1 main.cpp ${SOURCE} ${HEADERS})
target_include_directories(BecnhmarkCase1 PUBLIC ${common_dir})
# find_package(libInterpolate REQUIRED)
# find_package (Eigen3 3.3 REQUIRED)


target_link_libraries(BecnhmarkCase1 iMRS_library blas)
