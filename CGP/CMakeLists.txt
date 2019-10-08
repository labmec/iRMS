########### next target ###############

file(GLOB commonHeaders ../*.h)
file(GLOB commonSources ../*.cpp)

file(GLOB HEADERS *.h)
file(GLOB SOURCE *.cpp)

add_executable(CGP CGP.cpp ${lib_files} ${SOURCE} ${HEADERS} ${commonSources} ${commonHeaders})



target_include_directories(CGP PUBLIC ${lib_dir})

target_link_libraries(CGP pz ${OpenCV_LIBS})
file (COPY mathematica_notebooks/PlotGenerator.nb/ DESTINATION Debug/PlotGenerator.nb)