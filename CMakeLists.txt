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
find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../neopz_install/)

# Check if PZ was installed with MKL 
check_pz_opt(PZ_USING_MKL _PZ_HAS_MKL)
# TBD: mkl message not working. Need fixing on PZ
# message("mymessage = ${_PZ_HAS_MKL}")
# if(NOT _PZ_HAS_MKL)
#     message(FATAL_ERROR "PZ NOT INSTALED WITH MKL. Please do so!")
# endif()

# ----- Adding subdirectory of iMRS library -----
add_subdirectory(Common)

# ----- Subdirectories that are most used -----
add_subdirectory(BenchmarkCase1)
add_subdirectory(OnlyFractures)
add_subdirectory(2Elements2Spaces2D)

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