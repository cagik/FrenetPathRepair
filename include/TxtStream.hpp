#ifndef TXTSTREAM_HPP
#define TXTSTREAM_HPP
#include <fstream>
#include <vector>


template <typename T>
void DataToTxt(std::string txtPath, std::vector<std::pair<T, T>> &pts)
{
    std::fstream file(txtPath, std::ios::out);
    for(size_t i = 0; i < pts.size(); i++){
        file << pts[i].first << setprecision(10)  << "," << pts[i].second << setprecision(10) << std::endl;
    }
}

template <typename T>
void DataToTxt(std::string txtPath, std::vector<T> &pts)
{
    std::fstream file(txtPath, std::ios::out);
    for(size_t i = 0; i < pts.size(); i++){
        file << pts[i] << std::endl;
    }
}


#endif