########### next target ###############


file(GLOB common_headers *.h)
file(GLOB common_sources *.cpp)

set(common_files ${common_headers} ${common_sources} PARENT_SCOPE)
set(common_dir ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
