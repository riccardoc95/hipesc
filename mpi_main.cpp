//
// Created by Riccardo Ceccaroni on 31/01/25.
//

// mpic++ -std=c++17 -o my_program fastq.cpp paf.cpp jobs.cpp correction.cpp mpi_main.cpp
// mpirun -np 2 ./my_program

#include "Jobs.h"
#include "Correction.h"
#include <iostream>
#include <mpi.h>

#include <vector>
#include <string>

#include <chrono>


int main(int argc, char** argv) {
    //auto program_start_time = std::chrono::high_resolution_clock::now();


    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    Jobs jobs("data");
    auto job_list = jobs.get_job_list();
    size_t num_jobs = job_list.size();

    for (size_t i = world_rank; i < num_jobs; i += world_size) {
        std::cout << "Process " << world_rank << " processing job " << i + 1 << " of " << num_jobs << std::endl;
        auto job = jobs.get(job_list[i]);
        //print_job(job);
        std::cout << correction(job) << std::endl;
    }

    MPI_Finalize();


    // Record the end time for the entire program
    //auto program_end_time = std::chrono::high_resolution_clock::now();

    // Calculate the total duration (in seconds)
    //std::chrono::duration<double> program_duration = program_end_time - program_start_time;

    // Print the total elapsed time for the entire program
    //std::cout << "Total execution time: " << program_duration.count() << " seconds." << std::endl;

    return 0;
}