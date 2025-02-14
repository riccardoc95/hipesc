//
// Created by Riccardo Ceccaroni on 14/02/25.
//

#include "paf.h"


void parse_paf_line(std::string line_str, std::string& query_name, std::string& target_name, size_t& query_start, size_t& query_end, size_t& target_start, size_t& target_end) {
    const size_t pos1 = line_str.find('\t');
    const size_t pos2 = line_str.find('\t', pos1 + 1);
    const size_t pos3 = line_str.find('\t', pos2 + 1);
    const size_t pos4 = line_str.find('\t', pos3 + 1);
    const size_t pos5 = line_str.find('\t', pos4 + 1);
    const size_t pos6 = line_str.find('\t', pos5 + 1);
    const size_t pos7 = line_str.find('\t', pos6 + 1);
    const size_t pos8 = line_str.find('\t', pos7 + 1);
    const size_t pos9 = line_str.find('\t', pos8 + 1);

    query_name = line_str.substr(0, pos1);
    target_name = line_str.substr(pos5 + 1, pos6 - pos5 - 1);
    query_start = std::stoull(line_str.substr(pos2 + 1, pos3 - pos2 - 1));
    query_end = std::stoull(line_str.substr(pos3 + 1, pos4 - pos3 - 1));
    target_start = std::stoull(line_str.substr(pos7 + 1, pos8 - pos7 - 1));
    target_end = std::stoull(line_str.substr(pos8 + 1, pos9 - pos8 - 1));
}