#pragma once

#include <string>
#include <vector>

std::string trim(const std::string & source) {
    std::string s(source);
    s.erase(0,s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::vector<std::string> split_str(std::string line, char delimiter = ' '){
    std::vector<std::string> words;
    std::string tmp = "";

    int word_index = 0;
    for(char &c : line){
        if(c == delimiter){
            // next word
            word_index++;
            words.push_back(tmp);
            tmp = "";
        }
        else{
            tmp += c;
        }
    }

    words.push_back(tmp);

    return words;
}