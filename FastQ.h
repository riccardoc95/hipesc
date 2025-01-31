//
// Created by Riccardo Ceccaroni on 30/01/25.
//

#ifndef FASTQ_H
#define FASTQ_H
#include <unordered_map>
#include <iostream>

class FastQ {
public:
    explicit FastQ(const std::string& folder);
    std::unordered_map<std::string, std::string> get(const std::string& read_name);

//private:
    std::string folder;
    std::unordered_map<std::string, std::pair<std::string, std::streampos>> index;

    void index_fastq_files();
    void save_index();
    void load_index();
    bool is_fastq_file(const std::string& filename);
};

#endif // FASTQ_H
