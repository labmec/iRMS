########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable(MHM_Test main.cpp ${SOURCE} ${HEADERS} ${common_files})
target_include_directories(MHM_Test PUBLIC ${common_dir})
find_package(libInterpolate REQUIRED)

target_link_libraries(MHM_Test pz libInterpolate::Interpolate)
