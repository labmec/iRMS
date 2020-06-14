########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable(LinearTracer main.cpp ${SOURCE} ${HEADERS} ${common_files})
target_include_directories(LinearTracer PUBLIC ${common_dir})
find_package(libInterpolate REQUIRED)

target_link_libraries(LinearTracer pz libInterpolate::Interpolate)