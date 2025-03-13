//
// Created by Riccardo Ceccaroni on 12/02/25.
// clang++ --std=c++17 -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp -o openmp openmp.cpp -lz -lzstd -O2


#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>  // For std::remove
#include <sstream>    // For std::stringstream
#include <string>
#include <mpi.h>
#include <time.h>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include <zstd.h>  // Zstd header
#include <sys/resource.h>

#include "utils.h"
#include "consent.h"



using namespace std;


void print_memory_usage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "Peak memory usage: " << usage.ru_maxrss * 1e-06 << " MB" << std::endl;
}

// Structure to hold the start and end line for each job.

struct Target{
    string target_name;
    size_t query_start;  // The start line for the job
    size_t query_end;    // The end line for the job
    size_t target_start;  // The start line for the job
    size_t target_end;    // The end line for the job
    std::string strand;
};

struct Job {
    string query_name;
    vector<Target> jobs;
};

struct Alignment {
    int query_start;
    int query_end;
    std::string target_name;
    int target_start;
    int target_end;
};


std::string datetime() {
    // Get current time with milliseconds
    auto now = std::chrono::system_clock::now();
    auto now_time = std::chrono::system_clock::to_time_t(now);
    auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    // Convert to local time
    std::tm localTime = *std::localtime(&now_time);

    // Format time as string
    std::ostringstream timeStr;
    timeStr << std::put_time(&localTime, "%Y-%m-%d %H:%M:%S")
            << '.' << std::setfill('0') << std::setw(3) << now_ms.count();

    return timeStr.str();
}


void parse_line(string line_str, std::string& query_name, std::string& target_name, size_t& query_start, size_t& query_end, size_t& target_start, size_t& target_end, std::string& strand) {

    size_t pos1 = line_str.find('\t');
    size_t pos2 = line_str.find('\t', pos1 + 1);
    size_t pos3 = line_str.find('\t', pos2 + 1);
    size_t pos4 = line_str.find('\t', pos3 + 1);
    size_t pos5 = line_str.find('\t', pos4 + 1);
    size_t pos6 = line_str.find('\t', pos5 + 1);
    size_t pos7 = line_str.find('\t', pos6 + 1);
    size_t pos8 = line_str.find('\t', pos7 + 1);
    size_t pos9 = line_str.find('\t', pos8 + 1);

    query_name = line_str.substr(0, pos1);
    target_name = line_str.substr(pos5 + 1, pos6 - pos5 - 1);
    query_start = std::stoull(line_str.substr(pos2 + 1, pos3 - pos2 - 1));
    query_end = std::stoull(line_str.substr(pos3 + 1, pos4 - pos3 - 1));
    target_start = std::stoull(line_str.substr(pos7 + 1, pos8 - pos7 - 1));
    target_end = std::stoull(line_str.substr(pos8 + 1, pos9 - pos8 - 1));
    strand = "-";
}


// Function to compress a string using Zstd
string compress_zstd(const string& str) {
    size_t compressedSize = ZSTD_compressBound(str.size());
    vector<char> compressedData(compressedSize);

    // Compress the string
    compressedSize = ZSTD_compress(compressedData.data(), compressedData.size(), str.data(), str.size(), 3); // level 3 compression

    if (ZSTD_isError(compressedSize)) {
        throw runtime_error("Zstd Compression failed");
    }

    return string(compressedData.data(), compressedSize);
}

// Function to decompress a string using Zstd
string decompress_zstd(const string& compressedStr, size_t originalSize) {
    vector<char> decompressedData(originalSize);

    size_t decompressedSize = ZSTD_decompress(decompressedData.data(), decompressedData.size(), compressedStr.data(), compressedStr.size());

    if (ZSTD_isError(decompressedSize)) {
        throw runtime_error("Zstd Decompression failed");
    }

    return string(decompressedData.data(), decompressedSize);
}

// Function to trim newline characters from strings
void trim_newline(string &str) {
    str.erase(remove(str.begin(), str.end(), '\n'), str.end());
}

int main(int argc, char *argv[]) {

    std::string fastq_file = "data/test.fastq";
    std::string paf_file = "data/test.paf";

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
    bool done_reading = false;

    std::string d = datetime();
    std::cout << d << " - START" << ", Rank: " << rank << ", Size: " << size << endl;


    // Open the FASTQ file
    ifstream fastq(fastq_file);//input.fastq");
    if (!fastq.is_open()) {
        cerr << "Error opening FASTQ file!" << endl;
        return -1;
    }

    // Open the PAF file
    ifstream paf(paf_file);//input.fastq");
    if (!paf.is_open()) {
        cerr << "Error opening PAF file!" << endl;
        return -1;
    }

    // Create an unordered map to store FASTQ entries (header as key, sequence as value)
    unordered_map<string, pair<string, size_t>> fastq_dict;

    vector<Job> jobs;
    size_t paf_lines = 0;
    std::string line;
    while (std::getline(paf, line)) {
        paf_lines++;
    }
    paf.clear();
    paf.seekg(0, std::ios::beg);

    // LINES PER RANK
    size_t paf_lines_per_rank = (paf_lines / size) + 1;


    // Variables to store header, sequence, and quality scores
    string header, sequence, plus_line, quality_scores;
    string compressed_header, compressed_sequence;
    size_t sequence_length;

    int c_job = 0;

    // Reading the file and storing entries in the unordered map
    #pragma omp parallel
    {
        // Thread 0 will be responsible for reading and storing data in the map
        if (omp_get_thread_num() == 0) {
            d = datetime();
            std::cout << d << " - READ FASTQ, Rank: " << rank << std::endl;
            while (getline(fastq, header)) {
                trim_newline(header);
                header = header.substr(1);

                // Read the sequence line
                if (!getline(fastq, sequence)) break;
                trim_newline(sequence);

                // Read the plus line (not needed, so we just discard it)
                if (!getline(fastq, plus_line)) break;
                trim_newline(plus_line);

                // Read the quality scores line (not needed, so we discard it)
                if (!getline(fastq, quality_scores)) break;
                trim_newline(quality_scores);

                compressed_header = compress_zstd(header);
                compressed_sequence = compress_zstd(sequence);

                // Critical section to modify the shared unordered_map
                #pragma omp critical
                {
                    fastq_dict[compressed_header] = {compressed_sequence, sequence.size()};
                }
            }

        }
        // Synchronize threads to ensure reading is done before accessing the map
        #pragma omp barrier

        // Other threads can now read from the map (no modification, just reading)
        if (omp_get_thread_num() == 0) {
            d = datetime();
            std::cout << d << " - READ PAF & CREATE/PROCESS JOBS, Rank: " << rank << std::endl;

            size_t start_paf_line = rank * paf_lines_per_rank;
            size_t current_lines = 0;

            std::cout << "Rank: " << rank
                      << " Start paf line: " << start_paf_line
                      << ", Paf lines per rank: " << paf_lines_per_rank
                      << ", Number of lines: " << paf_lines
                      << std::endl;


            string query_name, target_name, back_query_name, strand;
            size_t query_start, query_end, target_start, target_end;
            bool end_job = false;

            while (current_lines < start_paf_line) {
                std::getline(paf, line);
                current_lines++;
            }
            current_lines = 0;

            vector<Target> targets;
            while (std::getline(paf, line)){
                current_lines++;
                if (current_lines >= paf_lines_per_rank && end_job) {
                    break;
                }

                parse_line(line, query_name, target_name, query_start, query_end, target_start, target_end, strand);
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
                    targets.push_back({compress_zstd(target_name), query_start, query_end, target_start, target_end, strand});
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
                    Job job = jobs.back();
                    //cout << "Thread: " << omp_get_thread_num() << " JOB: " << job.query_name << endl;
                    jobs.pop_back();
                    omp_unset_lock(&lock);

                    string query = job.query_name;
                    auto [fst, snd] = fastq_dict[query];
                    query = fst;
                    size_t query_length = snd;
                    query = decompress_zstd(query, query_length);

                    //cout << query << endl;


                    vector<Overlap> overlaps;
                    overlaps.clear();
                    // Iterate over targets within each job
                    for (const auto& target_list : job.jobs) {
                        string target = target_list.target_name;

                        auto [fst, snd] = fastq_dict[target];
                        target = fst;
                        size_t target_length = snd;
                        target = decompress_zstd(target, target_length);
                        //cout << target << " TARGET"<< endl;
                        //cout << query.size() << " " << target.size() << " " <<target_list.query_start << " " << target_list.query_end << " " << target_list.target_start << " " << target_list.target_end << endl;
                        overlaps.push_back({target, target_list.query_start, target_list.query_end, target_list.target_start, target_list.target_end, target_list.strand});
                    }

                    d = datetime();
                    std::cout << d << " - PRE_CORRECTION, Rank: " << rank << std::endl;
                    //consent_correction(query, overlaps);
                    d = datetime();
                    std::cout << d << " - POST_CORRECTION, Rank: " << rank << std::endl;


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
        d = datetime();
        std::cout << d << " - END, Rank: " << rank << std::endl;
        print_memory_usage();
    }

    omp_destroy_lock(&lock);

    // Finalize MPI
    MPI_Finalize();


    return 0;
}
