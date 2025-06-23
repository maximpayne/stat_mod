#include <vector>
#include <sstream>
#include <string>
#include <algorithm>

inline std::vector<double> parseNumbers(const std::string &line) {
    std::vector<double> nums;
    std::stringstream ss(line);
    std::string token;
    while(std::getline(ss, token, ',')) {
        std::stringstream ts(token);
        double val;
        if(ts >> val) nums.push_back(val);
    }
    return nums;
}
