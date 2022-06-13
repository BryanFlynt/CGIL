/**
 * \file       cgil2.cpp
 * \author     Bryan Flynt
 * \date       Jun 09, 2022
 * \copyright  Copyright (C) 2021 Bryan Flynt - All Rights Reserved
 */

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/Triangulation_hierarchy_2.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/natural_neighbor_coordinates_2.h>

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common.hpp"

// ============================================================================
//                                 Control
// ============================================================================

/// Maximum number of Variables at each Source Location
/**
 */
constexpr std::size_t MAXVAR = 8;

// ============================================================================
//                                 Types
// ============================================================================

using kernel                  = CGAL::Exact_predicates_inexact_constructions_kernel;
using vertex_base             = CGAL::Triangulation_vertex_base_2<kernel>;
using hierarchy_vertex_base   = CGAL::Triangulation_hierarchy_vertex_base_2<vertex_base>;
using face_base               = CGAL::Triangulation_face_base_2<kernel>;
using data_structure          = CGAL::Triangulation_data_structure_2<hierarchy_vertex_base, face_base>;
using triangulation           = CGAL::Delaunay_triangulation_2<kernel, data_structure>;
using hierarchy_triangulation = CGAL::Triangulation_hierarchy_2<triangulation>;
using point_type              = typename kernel::Point_2;
using coord_type              = typename kernel::FT;

using array_type     = Array<double, MAXVAR>;
using coord_map      = std::unordered_map<point_type, array_type>;
using value_accessor = CGAL::Data_access<coord_map>;

/// Read Source ASCII File
/**
 */
std::size_t read_source_file_2(const std::string& file_name, std::vector<point_type>& xyz,
                               std::vector<array_type>& var) {
    std::ifstream file(file_name);
    if (not file) {
        std::cerr << "ERROR: File Did No Open" << std::endl;
        std::cerr << "Filename: " << file_name << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Read
    std::stringstream file_stream;
    file_stream << file.rdbuf();
    file.close();
    
    // Parse Header
    std::size_t num_dim;
    std::size_t num_points;
    std::size_t num_variables;
    file_stream >> num_dim;
    if (2 != num_dim) {
        std::cerr << "ERROR: Wrong Number of Dimensions In File" << std::endl;
        std::cerr << "Expected 2 but got " << num_dim << std::endl;
        std::exit(EXIT_FAILURE);
    }
    file_stream >> num_points;
    file_stream >> num_variables;
    if (num_variables > MAXVAR) {
        std::cerr << "ERROR: Cannot Support More than " << MAXVAR << " Variables" << std::endl;
        std::cerr << "Requested " << num_variables << " Variables" << std::endl;
        std::exit(EXIT_FAILURE);
    }


    // Parse Data
    double x, y;
    xyz.resize(num_points);
    var.resize(num_points);
    for (std::size_t i = 0; i < num_points; ++i) {
        file_stream >> x;
        file_stream >> y;
        xyz[i] = point_type(x, y);
        for (std::size_t j = 0; j < num_variables; ++j) {
            file_stream >> var[i][j];
        }
    }
    return num_variables;
}

/// Read Target ASCII File
/**
 */
void read_target_file_2(const std::string& file_name, std::vector<point_type>& xyz) {
    std::ifstream file(file_name);
    if (not file) {
        std::cerr << "ERROR: File Did No Open" << std::endl;
        std::cerr << "Filename: " << file_name << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Read
    std::stringstream file_stream;
    file_stream << file.rdbuf();
    file.close();
    
    // Parse Header
    std::size_t num_dim;
    std::size_t num_points;
    file_stream >> num_dim;
    if (2 != num_dim) {
        std::cerr << "ERROR: Wrong Number of Dimensions In File" << std::endl;
        std::cerr << "Expected 2 but got " << num_dim << std::endl;
        std::exit(EXIT_FAILURE);
    }
    file_stream >> num_points;

    // Parse Data
    double x, y;
    xyz.resize(num_points);
    for (std::size_t i = 0; i < num_points; ++i) {
        file_stream >> x;
        file_stream >> y;
        xyz[i] = point_type(x, y);
    }
}

/// Write Target + Solution ASCII File
/**
 */
void write_target_file_2(const std::vector<point_type>& xyz, const std::size_t nvar,
                         const std::vector<array_type>& var) {
    std::cout << std::setw(10) << 2;
    std::cout << std::setw(10) << xyz.size();
    std::cout << std::setw(10) << nvar;
    std::cout << std::endl;

    for (std::size_t i = 0; i < xyz.size(); ++i) {
        for (std::size_t n = 0; n < 2; ++n) {
            std::cout << std::setw(15) << std::setprecision(8) << std::scientific << xyz[i][n] << " ";
        }
        for (std::size_t n = 0; n < nvar; ++n) {
            std::cout << std::setw(15) << std::setprecision(8) << std::scientific << var[i][n] << " ";
        }
        std::cout << std::endl;
    }
}

// ============================================================================
//                                  MAIN
// ============================================================================

int main(int argc, char* argv[]) {
    // Command line arguments
    std::string source_file_name = "source.txt";
    std::string target_file_name = "target.txt";
    if (argc > 1) {
        source_file_name = std::string(argv[1]);
    }
    if (argc > 2) {
        target_file_name = std::string(argv[2]);
    }

    // Reading Source Data
    std::vector<point_type> source_coords;
    std::vector<array_type> source_variables;
    const std::size_t nvar = read_source_file_2(source_file_name, source_coords, source_variables);

    // Reading Target Locations
    std::vector<point_type> target_coords;
    read_target_file_2(target_file_name, target_coords);

    // Build Delaunay Triangulation
    hierarchy_triangulation dt(source_coords.begin(), source_coords.end());

    // Build mapping from source_coords to source_variables
    coord_map point_to_value;
    for (auto i = 0; i < source_coords.size(); ++i) {
        point_to_value.insert(std::make_pair(source_coords[i], source_variables[i]));
    }

    // Search for Targets
    std::vector<array_type> target_variables(target_coords.size());
    for (std::size_t i = 0; i < target_coords.size(); ++i) {
        std::vector<std::pair<point_type, double>> coords;
        double norm = CGAL::natural_neighbor_coordinates_2(dt, target_coords[i], std::back_inserter(coords)).second;
        target_variables[i] =
            CGAL::linear_interpolation(coords.begin(), coords.end(), norm, value_accessor(point_to_value));
    }
    //     array_type sln;
    //     sln.fill(0);
    //     for(auto& coord : coords){
    //         sln += coord.second * point_to_value[coord.first];
    //     }
    //     target_variables[i] = sln / norm;
    // }

    // Write Solution to Screen
    write_target_file_2(target_coords, nvar, target_variables);

    return EXIT_SUCCESS;
}