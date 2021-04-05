########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable(GeoMeshAndPropsGenerator main.cpp ${SOURCE} ${HEADERS} ${common_files})
target_include_directories(GeoMeshAndPropsGenerator PUBLIC ${common_dir})
find_package(libInterpolate REQUIRED)
target_link_libraries(GeoMeshAndPropsGenerator NeoPZ::pz libInterpolate::Interpolate)
