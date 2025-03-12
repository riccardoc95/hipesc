//
// Created by Riccardo Ceccaroni on 14/02/25.
//

#ifndef FASTQ_H
#define FASTQ_H

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>


void read_fastq(std::ifstream& fastq, std::string& header, std::string& sequence, std::string& plus_line, std::string& quality_scores, std::string& compressed_header, std::string& compressed_sequence, std::unordered_map<std::string, std::pair<std::string, size_t>>& fastq_dict);

#endif //FASTQ_H
