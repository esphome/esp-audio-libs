cmake_minimum_required(VERSION 3.10)
project(flac_to_wav)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the esp-audio-libs as a subdirectory (going up two levels to the root)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../.. ${CMAKE_CURRENT_BINARY_DIR}/esp-audio-libs)

# Create the executable
add_executable(flac_to_wav src/flac_to_wav.cpp)

# Output the binary to the project root directory instead of build/
set_target_properties(flac_to_wav PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link with esp-audio-libs
target_link_libraries(flac_to_wav PRIVATE esp-audio-libs)

# Include directories
target_include_directories(flac_to_wav PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../include
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Add optimization flags
target_compile_options(flac_to_wav PRIVATE -O2)