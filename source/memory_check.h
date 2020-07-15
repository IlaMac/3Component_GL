#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

//#double process_memory_in_mb(const std::string &name){
//#    std::ifstream filestream("/proc/self/status");
//#    std::string line;
//#    while (std::getline(filestream, line)){
//#        std::istringstream is_line(line);
//#        std::string key;
//#        if (std::getline(is_line, key, ':')){
//#            if (key == name){
//#                std::string value_str;
//#                if (std::getline(is_line, value_str)) {
//#                    // Extract the number
//#                    std::string::size_type sz;   // alias of size_t
//#                    int value = std::stoi (value_str,&sz);
//#                    // Now we have the value in kb
//#                    return value/1024.0;
//#                }
//#            }
//#        }
//#    }
//#
//#    return -1.0;
//#}

double process_memory_in_mb(std::string_view name) {
    std::ifstream filestream("/proc/self/status");
    std::string   line;
    while(std::getline(filestream, line)) {
        std::istringstream is_line(line);
        std::string        key;
        if(std::getline(is_line, key, ':')) {
            if(key == name) {
                std::string value_str;
                if(std::getline(is_line, value_str)) {
                    // Filter non-digit characters
                    value_str.erase(std::remove_if(value_str.begin(), value_str.end(),
                                           []( auto const& c ) -> bool { return not std::isdigit(c); } ), value_str.end());
                    // Extract the number
                    long long value = 0;
                    try{
                        std::string::size_type sz; // alias of size_t
                        value = std::stoll(value_str, &sz);
                    }catch(const std::exception & ex){
			    printf("Could not read mem usage from /proc/self/status: Failed to parse string [%s]: %s \n", value_str.c_str(), ex.what());
                    }
                    // Now we have the value in kb
                    return static_cast<double>(value) / 1024.0;

                }
            }
        }
    }
    return -1.0;
}

