########### next target ###############

file(GLOB commonHeaders ../*.h)
file(GLOB commonSources ../*.cpp)

file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)

add_executable(LinearTracer main.cpp ${SOURCE} ${HEADERS} ${commonSources} ${commonHeaders})

target_link_libraries(LinearTracer pz)