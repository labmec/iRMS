file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)

add_executable(HdivCollapsedTest mainregular.cpp ${SOURCE} ${HEADERS})
target_link_libraries(HdivCollapsedTest pz)