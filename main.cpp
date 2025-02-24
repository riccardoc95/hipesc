//
// Created by Riccardo Ceccaroni on 12/02/25.
// mpic++ --std=c++17 -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp -o openmp openmp.cpp -lz -lzstd -O2
// mpic++ --std=c++17  -fopenmp -o openmp openmp.cpp -lz -lzstd -O2


#include <mpi.h>
#include <omp.h>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>

#include "src/fastq.h"
#include "src/paf.h"
#include "src/compression.h"
#include "src/utils.h"
#include "src/correction.h"


int main(int argc, char *argv[]) {
    std::string fastq_file = "../data/test.fastq";
    std::string paf_file = "../data/test.paf";

    // Check for command-line arguments and override the default files if provided
    if (argc > 1) {
        fastq_file = argv[1];  // First argument is the FASTQ file
    }
    if (argc > 2) {
        paf_file = argv[2];  // Second argument is the PAF file
    }

    // Display the file paths
    std::cout << "FASTQ file: " << fastq_file << std::endl;
    std::cout << "PAF file: " << paf_file << std::endl;

    // Initialize MPI
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get process rank
    MPI_Comm_size(MPI_COMM_WORLD, &size);  // Get number of processes


    // Initialize OMP
    omp_lock_t lock;
    omp_init_lock(&lock);

    print_time_message("START", rank);

    // Open files
    std::ifstream fastq(fastq_file);
    if (!fastq.is_open()) {
        std::cerr << "Error opening FASTQ file!" << std::endl;
        return -1;
    }
    std::ifstream paf(paf_file);
    if (!paf.is_open()) {
        std::cerr << "Error opening PAF file!" << std::endl;
        return -1;
    }

    // Count paf lines for split jobs on rank
    size_t paf_lines = 0;
    std::string line;
    while (std::getline(paf, line)) {
        paf_lines++;
    }
    paf.clear();
    paf.seekg(0, std::ios::beg);
    bool done_reading = false;

    // Paf lines per rank
    size_t paf_lines_per_rank = (paf_lines / size) + 1;

    #pragma omp parallel
    {
        std::unordered_map<std::string, std::pair<std::string, size_t>> fastq_dict;
        std::vector<Job> jobs;

        if (omp_get_thread_num() == 0) {
            std::string header;
            std::string sequence;
            std::string plus_line;
            std::string quality_scores;
            std::string compressed_sequence;
            std::string compressed_header;

            print_time_message("READ FASTQ", rank);
            read_fastq(fastq, header, sequence, plus_line, quality_scores, compressed_header, compressed_sequence, fastq_dict);
        }
        #pragma omp barrier

        if (omp_get_thread_num() == 0) {
            //TODO: SISTEMARE INIZIO LETTURA, SE RIGA PRECEDENTE E' UGUALE ALLA SUCCESSIVA, SKIP!!
            print_time_message("READ PAF", rank);

            size_t start_paf_line = rank * paf_lines_per_rank;
            size_t current_lines = 0;

            std::string query_name, target_name, back_query_name;
            size_t query_start, query_end, target_start, target_end;
            bool end_job = false;

            while (current_lines < start_paf_line) {
                std::getline(paf, line);
                current_lines++;
            }
            current_lines = 0;

            std::vector<Target> targets;
            while (std::getline(paf, line)){
                current_lines++;
                if (current_lines >= paf_lines_per_rank && end_job) {
                    break;
                }

                parse_paf_line(line, query_name, target_name, query_start, query_end, target_start, target_end);
                if (back_query_name.empty()) {
                  back_query_name = query_name;
                }
                if (query_name != back_query_name || line.empty()) {


                  omp_set_lock(&lock);
                  jobs.push_back({compress_zstd(back_query_name), targets});
                  omp_unset_lock(&lock);
                  back_query_name = query_name;

                  targets.clear();
                  end_job = true;
                } else{
                    targets.push_back({compress_zstd(target_name), query_start, query_end, target_start, target_end});
                    end_job = false;
                }

            }
            omp_set_lock(&lock);
            done_reading = true;
            omp_unset_lock(&lock);

        }
        else{
            while(true){
                omp_set_lock(&lock);
                if (!jobs.empty()) {
                    Job job = std::move(jobs.back());

                    jobs.pop_back();
                    omp_unset_lock(&lock);

                    std::string query = job.query_name;
                    auto [fst, snd] = fastq_dict[query];
                    query = fst;
                    size_t query_length = snd;
                    query = decompress_zstd(query, query_length);

                    std::vector<Overlap> targets;
                    targets.clear();
                    for (const auto&[target_name, query_start, query_end, target_start, target_end] : job.jobs) {
                        std::string target = target_name;

                        auto [fst, snd] = fastq_dict[target];
                        target = fst;
                        size_t target_length = snd;
                        target = decompress_zstd(target, target_length);
                        targets.push_back({
                            target.substr(target_start, target_end - target_start),
                            query_start,
                            query_end,
                            target_start,
                            target_end});


                    }
                    correction(query,targets, "default");

                } else if (done_reading) {
                    omp_unset_lock(&lock);
                    break;
                } else {
                    omp_unset_lock(&lock);
                }

            }
        }

    }
    #pragma omp barrier

    fastq.close();
    paf.close();

    if (omp_get_thread_num() == 0) {
        print_time_message("DONE", rank);
        print_memory_usage();
    }

    omp_destroy_lock(&lock);

    // Finalize MPI
    MPI_Finalize();


    return 0;
}
