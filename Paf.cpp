#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <vector>
#include <tuple>

#include "Paf.h"

namespace fs = std::filesystem;


Paf::Paf(const std::string& folder) : folder(folder) {
    if (fs::exists(folder + "/paf.index")) {
        load_index();
    } else {
        index_paf_file();
        save_index();
    }
}

void Paf::index_paf_file() {
    for (const auto& entry : fs::directory_iterator(folder)) {
        const auto& filename = entry.path().filename().string();
        if (is_paf_file(filename)) {
            std::ifstream paf(entry.path(), std::ios::in);
            if (paf.is_open()) {
                std::streampos position = paf.tellg();
                std::streampos position_prev = paf.tellg();
                std::string prev;
                while (true) {
                    std::string line;
                    if (!std::getline(paf, line)) break;  // End of file
                    std::istringstream stream(line);
                    std::string query_name;
                    stream >> query_name;

                    if (index.find(query_name) == index.end()) {
                        index[query_name] = {position};

                        if (!prev.empty() && prev != query_name) {
                            index[prev].push_back(position_prev);
                        }
                        prev = query_name;
                    }

                    position_prev = position;
                    position = paf.tellg();
                }
                index[prev].push_back(position_prev);
            }
        }
    }
}

void Paf::save_index() {
    std::ofstream idx(folder + "/paf.index");
    for (const auto& entry : index) {
        const auto& query_name = entry.first;
        const auto& positions = entry.second;
        for (size_t i = 0; i < positions.size(); i++) {
            idx << query_name << "\t" << positions[i] << "\t" << positions[i] << "\n";
        }
    }
}

void Paf::load_index() {
    std::ifstream idx(folder + "/paf.index");
    std::string line;
    while (std::getline(idx, line)) {
        std::istringstream stream(line);
        std::string query_name, start_str, end_str;
        std::getline(stream, query_name, '\t');
        std::getline(stream, start_str, '\t');
        std::getline(stream, end_str);

        index[query_name].push_back(std::stoi(start_str));
        index[query_name].push_back(std::stoi(end_str));
    }
}

std::unordered_map<std::string, std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string>>> Paf::get(const std::string& read_name) {
    auto it = index.find(read_name);
    std::unordered_map<std::string, std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string>>> sets;

    if (it != index.end()) {
        // iterate over positions for the query
        const auto& positions = it->second;
        for (size_t i = 0; i < positions.size(); i += 2) {
            // Start position and end position
            std::streampos start = positions[i];
            std::streampos end = positions[i + 1];

            for (const auto& entry : fs::directory_iterator(folder)) {
                const auto& filename = entry.path().filename().string();
                if (is_paf_file(filename)) {
                    std::ifstream paf(entry.path(), std::ios::in);
                    if (paf.is_open()) {
                        paf.seekg(start);
                        while (paf.tellg() < end) {
                            std::string line;
                            if (!std::getline(paf, line)) break;
                            std::istringstream stream(line);
                            std::string query_name, query_start, query_end, target_name, target_start, target_end, _;
                            stream >> query_name >> _ >>query_start >> query_end >> _ >>target_name >> _ >> target_start >> target_end;

                            sets[query_name].emplace_back(target_name, query_start, query_end, target_start, target_end);
                        }
                    }
                }
            }
        }
    }
    return sets;
}


bool Paf::is_paf_file(const std::string& filename) {
    size_t len = filename.length();
    return (len >= 4 && filename.substr(len - 4) == ".paf");
}



