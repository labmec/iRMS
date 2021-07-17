########### SimpleCaseFracture ###############

file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)

add_executable( SimpleCaseFracture main.cpp ${SOURCE} ${HEADERS})
target_link_libraries( SimpleCaseFracture PUBLIC iMRS_library)
