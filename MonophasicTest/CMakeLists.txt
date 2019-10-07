########### next target ###############

file(GLOB HEADERS *.h)

file(GLOB SOURCE *.cpp)

add_executable(MonophasicTest main.cpp ${SOURCE} ${HEADERS})

target_link_libraries(MonophasicTest pz)