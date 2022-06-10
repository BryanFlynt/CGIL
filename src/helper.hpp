/**
 * \file       helper.hpp
 * \author     Bryan Flynt
 * \date       Jun 08, 2022
 * \copyright  Copyright (C) 2021 Bryan Flynt - All Rights Reserved
 */
#pragma once

#include <array>

#include <vector>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>

#include <cstdint>

template<typename T>
inline void 
read_sources(const std::string& file_name, 
             std::vector<std::pair<double,double>>& source_coordinates, 
             std::vector<T>& source_variables) {

    std::ifstream file(file_name);
    if(not file ){
        std::cerr << "ERROR: File Did No Open" << std::endl;
        std::cerr << "Filename: " << file_name << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Get size of File
    file.seekg(0,std::ios::end);
    std::streampos length = file.tellg();
    file.seekg(0,std::ios::beg);

    // Allocate std::vector as Buffer and read
    std::vector<char> buffer(length);
    file.read(buffer.data(),length);

    // Create string stream
    std::stringstream localStream;
    localStream.rdbuf()->pubsetbuf(buffer.data(),length);

    // Parse Header
    std::size_t num_points;
    std::size_t num_variables;
    localStream >> num_points;
    localStream >> num_variables;

    // Parse Data
    double x, y, v;
    source_coordinates.resize(num_points);
    source_variables.resize(num_points);
    if(source_variables[0].size() < num_variables){
        std::cerr << "Need to increase the size of Variable Array" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    for(auto i = 0; i < num_points; ++i){
        localStream >> x;
        localStream >> y;
        source_coordinates[i] = std::make_pair(x,y);  
        for(auto j = 0; j < num_variables; ++j){
            localStream >> v;
            source_variables[i][j] = v;
        }
    }
}



inline void 
read_targets(const std::string& file_name, 
             std::vector<std::pair<double,double>>& target_coordinates) {

    std::ifstream file(file_name);
    if(not file ){
        std::cerr << "ERROR: File Did No Open" << std::endl;
        std::cerr << "Filename: " << file_name << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Get size of File
    file.seekg(0,std::ios::end);
    std::streampos length = file.tellg();
    file.seekg(0,std::ios::beg);

    // Allocate std::vector as Buffer and read
    std::vector<char> buffer(length);
    file.read(buffer.data(),length);

    // Create string stream
    std::stringstream localStream;
    localStream.rdbuf()->pubsetbuf(buffer.data(),length);

    // Parse Header
    std::size_t num_points;
    localStream >> num_points;

    // Parse Data
    double x, y;
    target_coordinates.resize(num_points);
    for(auto i = 0; i < num_points; ++i){
        localStream >> x;
        localStream >> y;
        target_coordinates[i] = std::make_pair(x,y);  
    }
}


template<typename T, std::size_t N>
std::array<T,N> operator+(const std::array<T,N>& a, const std::array<T,N>& b){
    std::array<T,N> ans(a);
    for(std::size_t i = 0; i < a.size(); ++i){
        ans[i] += b[i];
    }    
    return ans;
}

template<typename T, std::size_t N>
std::array<T,N> operator-(const std::array<T,N>& a, const std::array<T,N>& b){
    std::array<T,N> ans(a);
    for(std::size_t i = 0; i < a.size(); ++i){
        ans[i] -= b[i];
    }    
    return ans;
}

template<typename T, std::size_t N>
std::array<T,N> operator*(const std::array<T,N>& a, const std::array<T,N>& b){
    std::array<T,N> ans(a);
    for(std::size_t i = 0; i < a.size(); ++i){
        ans[i] *= b[i];
    }    
    return ans;
}

template<typename T, std::size_t N>
std::array<T,N> operator/(const std::array<T,N>& a, const std::array<T,N>& b){
    std::array<T,N> ans(a);
    for(std::size_t i = 0; i < a.size(); ++i){
        ans[i] /= b[i];
    }    
    return ans;
}

template<typename T, std::size_t N>
std::array<T,N> operator*(const T a, const std::array<T,N>& b){
    std::array<T,N> ans(b);
    for(std::size_t i = 0; i < a.size(); ++i){
        ans[i] *= a;
    }    
    return ans;
}

template<typename T, std::size_t N>
std::array<T,N> operator*(const std::array<T,N>& a, const T b){
    return b * a;
}


template<typename T, std::size_t N>
std::array<T,N> operator/(const T a, const std::array<T,N>& b){
    std::array<T,N> ans;
    for(std::size_t i = 0; i < a.size(); ++i){
        ans[i] = a / b[i];
    }    
    return ans;
}

template<typename T, std::size_t N>
std::array<T,N> operator/(const std::array<T,N>& a, const T b){
    std::array<T,N> ans;
    for(std::size_t i = 0; i < a.size(); ++i){
        ans[i] = a[i] / b;
    }    
    return ans;
}