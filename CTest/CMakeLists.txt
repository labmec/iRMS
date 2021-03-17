# Code testing. Please run these before Pushing or PR
enable_testing()

# SimpleCaseFractureTwoElements
add_test(NAME SimpleCaseFractureTwoElements 
        COMMAND SimpleCaseFractureTwoElements 
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/SimpleCaseFractureTwoElements/Debug"    # Jose uses this /Debug directory
        # WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/SimpleCaseFractureTwoElements"    # But I guess this would be the default
)