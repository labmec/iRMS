########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable( OnlyTwoFractures main.cpp ${SOURCE} ${HEADERS})
target_include_directories( OnlyTwoFractures PUBLIC ${common_dir})
# find_package(libInterpolate REQUIRED)
# find_package (Eigen3 3.3 REQUIRED)

find_package(BLAS REQUIRED)
target_link_libraries( OnlyTwoFractures iMRS_library ${BLAS_LIBRARIES})
