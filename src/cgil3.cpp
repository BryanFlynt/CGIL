/**
 * \file       cgil3.cpp
 * \author     Bryan Flynt
 * \date       Jun 09, 2022
 * \copyright  Copyright (C) 2021 Bryan Flynt - All Rights Reserved
 */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/squared_distance_3.h> //for 3D functions
#include <CGAL/interpolation_functions.h>


#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <limits>

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
using vertex_base             = CGAL::Triangulation_vertex_base_3<kernel>;
using cell_base               = CGAL::Delaunay_triangulation_cell_base_3<kernel>;
using data_structure          = CGAL::Triangulation_data_structure_3<vertex_base, cell_base>;
using triangulation           = CGAL::Delaunay_triangulation_3<kernel, data_structure, CGAL::Fast_location>;
using point_type              = typename triangulation::Point;
using cell_handle             = typename triangulation::Cell_handle;
using vertex_handle           = typename triangulation::Vertex_handle;

using array_type     = Array<double, MAXVAR>;
using coord_map      = std::unordered_map<point_type, array_type>;
using value_accessor = CGAL::Data_access<coord_map>;


/// Read Source ASCII File
/**
 */
std::size_t read_source_file_3(const std::string& file_name, std::vector<point_type>& xyz,
                               std::vector<array_type>& var) {
    std::ifstream file(file_name);
    if (not file) {
        std::cerr << "ERROR: File Did No Open" << std::endl;
        std::cerr << "Filename: " << file_name << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Get size of File
    file.seekg(0, std::ios::end);
    std::streampos length = file.tellg();
    file.seekg(0, std::ios::beg);

    // Allocate std::vector as Buffer and read
    std::vector<char> buffer(length);
    file.read(buffer.data(), length);

    // Create string stream
    std::stringstream localStream;
    localStream.rdbuf()->pubsetbuf(buffer.data(), length);

    // Parse Header
    std::size_t num_dim;
    std::size_t num_points;
    std::size_t num_variables;
    localStream >> num_dim;
    if (3 != num_dim) {
        std::cerr << "ERROR: Wrong Number of Dimensions In File" << std::endl;
        std::cerr << "Expected 3 but got " << num_dim << std::endl;
        std::exit(EXIT_FAILURE);
    }
    localStream >> num_points;
    localStream >> num_variables;
    if (num_variables > MAXVAR) {
        std::cerr << "ERROR: Cannot Support More than " << MAXVAR << " Variables" << std::endl;
        std::cerr << "Requested " << num_variables << " Variables" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Parse Data
    double x, y, z;
    xyz.resize(num_points);
    var.resize(num_points);
    for (std::size_t i = 0; i < num_points; ++i) {
        localStream >> x;
        localStream >> y;
        localStream >> z;
        xyz[i] = point_type(x, y, z);
        for (std::size_t j = 0; j < num_variables; ++j) {
            localStream >> var[i][j];
        }
    }
    return num_variables;
}

/// Read Target ASCII File
/**
 */
void read_target_file_3(const std::string& file_name, std::vector<point_type>& xyz) {
    std::ifstream file(file_name);
    if (not file) {
        std::cerr << "ERROR: File Did No Open" << std::endl;
        std::cerr << "Filename: " << file_name << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Get size of File
    file.seekg(0, std::ios::end);
    std::streampos length = file.tellg();
    file.seekg(0, std::ios::beg);

    // Allocate std::vector as Buffer and read
    std::vector<char> buffer(length);
    file.read(buffer.data(), length);

    // Create string stream
    std::stringstream localStream;
    localStream.rdbuf()->pubsetbuf(buffer.data(), length);

    // Parse Header
    std::size_t num_dim;
    std::size_t num_points;
    localStream >> num_dim;
    if (3 != num_dim) {
        std::cerr << "ERROR: Wrong Number of Dimensions In File" << std::endl;
        std::cerr << "Expected 3 but got " << num_dim << std::endl;
        std::exit(EXIT_FAILURE);
    }
    localStream >> num_points;

    // Parse Data
    double x, y, z;
    xyz.resize(num_points);
    for (std::size_t i = 0; i < num_points; ++i) {
        localStream >> x;
        localStream >> y;
        localStream >> z;
        xyz[i] = point_type(x, y, z);
    }
}

/// Write Target + Solution ASCII File
/**
 */
void write_target_file_3(const std::vector<point_type>& xyz, const std::size_t nvar,
                         const std::vector<array_type>& var) {
    std::cout << std::setw(10) << 3;
    std::cout << std::setw(10) << xyz.size();
    std::cout << std::setw(10) << nvar;
    std::cout << std::endl;

    for (std::size_t i = 0; i < xyz.size(); ++i) {
        for (std::size_t n = 0; n < 3; ++n) {
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
    const std::size_t nvar = read_source_file_3(source_file_name, source_coords, source_variables);

    // Reading Target Locations
    std::vector<point_type> target_coords;
    read_target_file_3(target_file_name, target_coords);

    // Build Delaunay Triangulation
    triangulation dt(source_coords.begin(), source_coords.end());

    // Build mapping from source_coords to source_variables
    coord_map point_to_value;
    for (auto i = 0; i < source_coords.size(); ++i) {
        point_to_value.insert(std::make_pair(source_coords[i], source_variables[i]));
    }


    // Search for Targets
    cell_handle c;
    std::vector<array_type> target_variables(target_coords.size());
    for (std::size_t i = 0; i < target_coords.size(); ++i) {

        // Assuming some structure to the target points the 
        // previous cell might be a good first guess ???
        c = dt.locate(target_coords[i], c);

        // Methods do not exists within CGAL for 3D Weighting
        // We use a basic Inverse Distance Weighting functions
        // - Find distance^2 for each vertex to target
        // - If Exact Match then Done
        // - Else then calculate weights
        std::size_t n;
        std::array<double,4> sqdist;
        for(n = 0; n < 4; ++n) {
            sqdist[n] = CGAL::squared_distance(target_coords[i], c->vertex(n)->point());
            if( sqdist[n] < std::numeric_limits<double>::epsilon() ){
                break;
            }
        }

        if( n < 4 ){
            target_variables[i] = point_to_value[c->vertex(n)->point()];
        }
        else {
            double norm = 0;
            std::vector<std::pair<point_type, double>> coords(4);
            for(std::size_t n = 0; n < 4; ++n) {
                auto w = static_cast<double>(1) / sqdist[n];
                norm += w;
                coords[n] = std::make_pair(c->vertex(n)->point(), w);
            }
            target_variables[i] = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, value_accessor(point_to_value));
        }
    }

    // // Write Solution to Screen
    write_target_file_3(target_coords, nvar, target_variables);

    return EXIT_SUCCESS;
}