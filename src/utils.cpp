//
// Created by Riccardo Ceccaroni on 14/02/25.
//

#include "utils.h"

#include <iostream>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>


void print_memory_usage() {
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    std::cout << "Peak memory usage: " << static_cast<double>(usage.ru_maxrss) * 1e-06 << " MB" << std::endl;
}

std::string datetime() {
    const auto now = std::chrono::system_clock::now();
    const auto  now_time = std::chrono::system_clock::to_time_t(now);
    const auto  now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

    const std::tm localTime = *std::localtime(&now_time);
    std::ostringstream timeStr;
    timeStr << std::put_time(&localTime, "%Y-%m-%d %H:%M:%S")
            << '.' << std::setfill('0') << std::setw(3) << now_ms.count();

    return timeStr.str();
}

void trim_newline(std::string &str) {
    str.erase(std::remove(str.begin(), str.end(), '\n'), str.end());
}

void print_time_message(const std::string& message, const int rank) {
    const std::string d = datetime();
    std::cout << d << " - " << message << ", Rank: " << rank << std::endl;
}


