cmake_minimum_required(VERSION 3.29)
project(hipes)

# Enable MPI
find_package(MPI REQUIRED)

# Include directories
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 17)

add_executable(hipes main.cpp
        fastq.cpp
        fastq.h
        paf.cpp
        paf.h
        jobs.cpp
        jobs.h
        mpi_main.cpp
        Correction.cpp
        Correction.h)

# Link MPI libraries
target_link_libraries(hipes ${MPI_LIBRARIES})
