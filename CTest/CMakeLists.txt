# Code testing. Please run these before Pushing or PR
enable_testing()

# MHMSimple
add_test(NAME MHMSimple COMMAND MHMSimple WORKING_DIRECTORY ${CMAKE_BUILD_DIR}/MHMSimple/Debug)