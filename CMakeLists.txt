cmake_minimum_required(VERSION 3.10)
project(hipesc LANGUAGES CXX)

find_package(MPI REQUIRED)
find_package(Threads REQUIRED)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Find OpenMP
find_package(OpenMP)
if(NOT OpenMP_CXX_FOUND)
    message(WARNING "OpenMP not found via CMake, manually setting paths...")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_CXX_LIB_NAMES "omp")
endif()

# zstd
execute_process(COMMAND brew --prefix zstd OUTPUT_VARIABLE ZSTD_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
include_directories(${ZSTD_PREFIX}/include)
link_directories(${ZSTD_PREFIX}/lib)
#target_link_libraries(hipesc PRIVATE zstd)

# Apply OpenMP flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L$(brew --prefix libomp)/lib -lomp")

# Get libomp paths from Homebrew
execute_process(COMMAND brew --prefix libomp OUTPUT_VARIABLE LIBOMP_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
include_directories(${LIBOMP_PREFIX}/include)
link_directories(${LIBOMP_PREFIX}/lib)

# Add executable
add_executable(hipesc main.cpp
        src/utils.cpp
        src/utils.h
        src/compression.cpp
        src/compression.h
        src/correction.cpp
        src/correction.h
        src/paf.cpp
        src/paf.h
        src/fastq.cpp
        src/fastq.h
)

# Link required libraries
target_link_libraries(hipesc PRIVATE MPI::MPI_CXX z zstd)

# Set compile options
target_compile_options(hipesc PRIVATE -O2)
