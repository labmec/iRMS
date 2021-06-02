########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable( SimpleCaseTwoFractures main.cpp ${SOURCE} ${HEADERS})
target_include_directories( SimpleCaseTwoFractures PUBLIC ${common_dir})
# find_package(libInterpolate REQUIRED)
# find_package (Eigen3 3.3 REQUIRED)

find_package(BLAS REQUIRED)
target_link_libraries( SimpleCaseTwoFractures iMRS_library ${BLAS_LIBRARIES})
