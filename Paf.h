//
// Created by Riccardo Ceccaroni on 30/01/25.
//

#ifndef PAF_H
#define PAF_H
#include <unordered_map>
#include <iostream>

class Paf {
public:
    explicit Paf(const std::string& folder);
    std::unordered_map<std::string, std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string>>> get(const std::string& read_name);

//private:
    std::string folder;
    std::unordered_map<std::string, std::vector<std::streampos>> index;

    void index_paf_file();
    void save_index();
    void load_index();
    bool is_paf_file(const std::string& filename);
};

#endif // PAF_H
