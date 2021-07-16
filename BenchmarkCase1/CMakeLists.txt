########### BenchmarkCase1 ###############

file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)

add_executable(BenchmarkCase1 main.cpp ${SOURCE} ${HEADERS})
target_link_libraries(BenchmarkCase1 PUBLIC iMRS_library)
