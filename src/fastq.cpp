//
// Created by Riccardo Ceccaroni on 14/02/25.
//

#include "fastq.h"

#include "compression.h"
#include "utils.h"


void read_fastq(std::ifstream& fastq, std::string& header, std::string& sequence, std::string& plus_line, std::string& quality_scores, std::string& compressed_header, std::string& compressed_sequence, std::unordered_map<std::string, std::pair<std::string, size_t>>& fastq_dict) {
    while (std::getline(fastq, header)) {
        trim_newline(header);
        header = header.substr(1);

        if (!getline(fastq, sequence)) break;
        trim_newline(sequence);

        if (!getline(fastq, plus_line)) break;
        trim_newline(plus_line);

        if (!getline(fastq, quality_scores)) break;
        trim_newline(quality_scores);

        compressed_header = compress_zstd(header);
        compressed_sequence = compress_zstd(sequence);

        #pragma omp critical
        {
            fastq_dict[compressed_header] = {compressed_sequence, sequence.size()};
        }
    }
}