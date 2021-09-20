########### next target ###############

# file(GLOB HEADERS *.h)
# file(GLOB SOURCE *.cpp)

# add_executable( SimpleCaseTwoFractures main.cpp ${SOURCE} ${HEADERS})
add_executable( SimpleCaseTwoFractures main.cpp)
target_link_libraries( SimpleCaseTwoFractures PUBLIC iMRS_library)
