########### next target ###############

file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)

add_executable(LinearTracer main.cpp ${SOURCE} ${HEADERS} ${common_files})
target_include_directories(LinearTracer PUBLIC ${common_dir})

target_link_libraries(LinearTracer pz)
