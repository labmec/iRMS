cmake_minimum_required(VERSION 3.14)

# ----- Basic properties -----
project (iMRS LANGUAGES CXX)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ----- Linking with accelerate lib if on macos -----
IF(APPLE)
   add_definitions(-DMACOSX)
   set(MACOSX MACOSX)
   FIND_LIBRARY(ACCELERATE_LIB Accelerate)
   link_libraries(${ACCELERATE_LIB})
ENDIF(APPLE)

# ----- Find PZ. By default searches one directory up -----
find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/ ${CMAKE_SOURCE_DIR}/neopz_install/)
 

# ----- Check if PZ was installed with MKL -----
check_pz_opt(PZ_USING_MKL _PZ_HAS_MKL)
if(${_PZ_HAS_MKL} STREQUAL "OFF")
    message(FATAL_ERROR "\nERROR! NEOPZ WITHOUT MKL IS NOT SUPPORTED. PLEASE COMPILE NEOPZ WITH MKL")
endif()


file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
)
add_custom_target( GlobalFiles SOURCES ${PZ_SRC})

# ----- To enable unit test using Catch2 UnitTest framework -----
option(BUILD_UNITTESTING "Whether the Catch2 UnitTest will be linked in" OFF)
if(BUILD_UNITTESTING)
  include(cmake/add_unit_test.cmake)
  include(cmake/EnableCatch2.cmake)
  enable_catch2()
endif()


# ----- Adding subdirectory of DFN library -----
# linking with GMSH
find_library(GMSH_LIB gmsh HINTS "/home/pedro/projetos/library/gmsh-install/lib" "../lib/gmsh-install/lib" "../gmsh-4.8.4-MacOSX-sdk/lib" NO_DEFAULT_PATH)
find_path(GMSH_INC gmsh.h HINTS "/home/pedro/projetos/library/gmsh-install/include" "../lib/gmsh-install/include" "../gmsh-4.8.4-MacOSX-sdk/include" NO_DEFAULT_PATH)

link_libraries(${GMSH_LIB})
if(NOT GMSH_LIB)
  message(FATAL_ERROR "Could not find libgmsh")
endif()

if(NOT GMSH_INC)
  message(FATAL_ERROR "Could not find gmsh.h")
endif()
add_subdirectory(src_dfn)

# ----- Adding subdirectory of iMRS library -----
add_subdirectory(Common)


# ----- Subdirectories that are most used -----
add_subdirectory(BenchmarkCase1)
add_subdirectory(OnlyFractures)
add_subdirectory(2Elements2Spaces2D)
# add_subdirectory(SimpleCaseTwoFractures)
add_subdirectory(SimpleCaseOneFracture)
add_subdirectory(CubicDomainNoFractures)
add_subdirectory(DFNIMRS)


# ----- Enables testing -----
if(BUILD_UNITTESTING)
    add_subdirectory(UnitTests)    
    ENABLE_TESTING()
    include(CTest)
endif()

# ----- Subdirectories that need fixing to new pz (with the crazy materials) -----
# add_subdirectory(SimpleCaseFracture)
# add_subdirectory(SimpleCaseTwoFractures)
# add_subdirectory(SimpleCaseThreeFractures)
# add_subdirectory(HDivSimple)

# ----- Subdirectories that need fixing and are seldom used -----
# add_subdirectory(MHMSimple)
# add_subdirectory(GeoMeshAndPropsGenerator)
# add_subdirectory(CGP)
# add_subdirectory(HdivCollapsed)
# add_subdirectory(Paper2D)