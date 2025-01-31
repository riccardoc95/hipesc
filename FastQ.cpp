//
// Created by Riccardo Ceccaroni on 30/01/25.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <filesystem>

#include "FastQ.h"


namespace fs = std::filesystem;

FastQ::FastQ(const std::string& folder) : folder(folder) {
    if (fs::exists(folder + "/fastq.index")) {
        load_index();
    } else {
        index_fastq_files();
        save_index();
    }
}

void FastQ::index_fastq_files() {
    for (const auto& entry : fs::directory_iterator(folder)) {
        const auto& filename = entry.path().filename().string();
        if (is_fastq_file(filename)) {
            std::ifstream fq(entry.path(), std::ios::in);
            if (fq.is_open()) {
                std::streampos position = fq.tellg();
                std::string header;
                while (std::getline(fq, header)) {
                    if (header.empty()) break;

                    std::string name = header.substr(1, header.find(" ") - 1);  // Remove '@' and extract name
                    index[name] = {filename, position};

                    // Skip the next three lines (sequence, '+', quality)
                    std::getline(fq, header);
                    std::getline(fq, header);
                    std::getline(fq, header);

                    position = fq.tellg();
                }
            }
        }
    }
}

void FastQ::save_index() {
    std::ofstream idx(folder + "/fastq.index");
    for (const auto& entry : index) {
        const auto& name = entry.first;
        const auto& [filename, position] = entry.second;
        idx << name << "\t" << filename << "\t" << position << "\n";
    }
}

void FastQ::load_index() {
    std::ifstream idx(folder + "/fastq.index");
    std::string line;
    while (std::getline(idx, line)) {
        std::istringstream stream(line);
        std::string name, filename;
        std::string position_str;
        std::getline(stream, name, '\t');
        std::getline(stream, filename, '\t');
        std::getline(stream, position_str);

        index[name] = {filename, std::stoi(position_str)};
    }
}

std::unordered_map<std::string, std::string> FastQ::get(const std::string& read_name) {
    auto it = index.find(read_name);
    if (it != index.end()) {
        const auto& [filename, position] = it->second;
        std::ifstream fq(folder + "/" + filename, std::ios::in);
        if (fq.is_open()) {
            fq.seekg(position);
            std::string header, sequence, plus, quality;
            std::getline(fq, header);
            std::getline(fq, sequence);
            std::getline(fq, plus);
            std::getline(fq, quality);

            std::unordered_map<std::string, std::string> result;
            result["read_name"] = header.substr(1);  // Remove '@'
            result["sequence"] = sequence;
            result["plus"] = plus;
            result["quality"] = quality;

            return result;
        }
    }
    return {};
}

bool FastQ::is_fastq_file(const std::string& filename) {
    size_t len = filename.length();
    return (len >= 5 && filename.substr(len - 5) == ".fastq") ||
           (len >= 3 && filename.substr(len - 3) == ".fq");
}
