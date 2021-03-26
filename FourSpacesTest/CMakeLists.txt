########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable(FourSpacesTest main.cpp ${SOURCE} ${HEADERS} ${common_files})
target_include_directories(FourSpacesTest PUBLIC ${common_dir})
find_package(libInterpolate REQUIRED)

target_link_libraries(FourSpacesTest NeoPZ::pz libInterpolate::Interpolate)