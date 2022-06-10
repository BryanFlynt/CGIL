/**
 * \file       make_random_sources.cpp
 * \author     Bryan Flynt
 * \date       Jun 09, 2022
 * \copyright  Copyright (C) 2021 Bryan Flynt - All Rights Reserved
 */

#include "common.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

int main(int argc, char *argv[]) {

    std::size_t num_dimensions = 2;
    std::size_t num_points     = 10;
    std::size_t num_variables  = 2;

    if( argc > 1 ){
        num_dimensions = std::atoi(argv[1]);
    }
    if( argc > 2 ){
        num_points = std::atoi(argv[2]);
    }
    if( argc > 3 ){
        num_variables = std::atoi(argv[3]);
    }

    const std::size_t min_points = std::pow(2,num_dimensions) + 1;
    if( num_points < min_points ){
        std::cerr << "ERROR: Need at least " << min_points << "points" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Write Header
    std::cout << std::setw(10) << num_dimensions;
    std::cout << std::setw(10) << num_points;
    std::cout << std::setw(10) << num_variables;
    std::cout << std::endl;

    // Set Bounds of Region
    const double xmin    = -1000;
    const double xmax    = +1000;
    const auto xyz_range = std::make_pair(xmin,xmax);

    // Allocate Buffers
    std::vector<double> xyz(num_dimensions);
    std::vector<double> var(num_variables);

    // Write Corners
    if(num_dimensions == 2){
        xyz = { 0.0,  0.0}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmin, xmin}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmin, xmax}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmax, xmax}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmin, xmax}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
    }
    else if(num_dimensions == 3){
        xyz = { 0.0,  0.0,  0.0}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmin, xmin, xmin}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmin, xmax, xmin}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmax, xmax, xmin}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmin, xmax, xmin}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmin, xmin, xmax}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmin, xmax, xmax}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmax, xmax, xmax}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
        xyz = {xmin, xmax, xmax}; common::fill_variables(xyz,var); common::dump_line(xyz, var);
    }
    else {
        std::cerr << "ERROR: Only Support 2 & 3 Dimensions" << std::endl;
        exit(EXIT_FAILURE);
    }

    for(std::size_t i = min_points; i < num_points; ++i){
        common::fill_random(xyz_range, xyz);
        common::fill_variables(xyz,var);
        common::dump_line(xyz, var);
    }

    return 0;
}