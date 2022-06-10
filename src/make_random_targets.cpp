/**
 * \file       make_random_targets.cpp
 * \author     Bryan Flynt
 * \date       Jun 09, 2022
 * \copyright  Copyright (C) 2021 Bryan Flynt - All Rights Reserved
 */

#include <cstdlib>
#include <iomanip>
#include <iostream>

#include "common.hpp"

int main(int argc, char *argv[]) {
    std::size_t num_dimensions = 2;
    std::size_t num_points     = 10;

    if (argc > 1) {
        num_dimensions = std::atoi(argv[1]);
    }
    if (argc > 2) {
        num_points = std::atoi(argv[2]);
    }

    const std::size_t min_points = 1;
    if (num_points < min_points) {
        std::cerr << "ERROR: Need at least " << min_points << "points" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Write Header
    std::cout << std::setw(10) << num_dimensions;
    std::cout << std::setw(10) << num_points;
    std::cout << std::endl;

    // Set Bounds of Region
    const double xmin    = -10;
    const double xmax    = +10;
    const auto xyz_range = std::make_pair(xmin, xmax);

    // Write Center for Easy Check
    std::vector<double> xyz(num_dimensions, 0.0);
    common::dump_line(xyz);

    // Write Random Points
    for (std::size_t i = 1; i < num_points; ++i) {
        common::fill_random(xyz_range, xyz);
        common::dump_line(xyz);
    }

    return 0;
}