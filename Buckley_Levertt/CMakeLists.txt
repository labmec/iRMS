########### next target ###############



file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)



add_executable(Buckley_Levertt main.cpp ${SOURCE} ${HEADERS} ${common_files})
target_include_directories(Buckley_Levertt PUBLIC ${common_dir})

target_link_libraries(Buckley_Levertt pz)