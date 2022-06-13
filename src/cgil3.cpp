/**
 * \file       cgil3.cpp
 * \author     Bryan Flynt
 * \date       Jun 09, 2022
 * \copyright  Copyright (C) 2021 Bryan Flynt - All Rights Reserved
 */

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/interpolation_functions.h>
#include <CGAL/squared_distance_3.h>  //for 3D functions

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common.hpp"

inline constexpr bool PRE_ALLOCATE_BUFFER = true;

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

using kernel         = CGAL::Exact_predicates_inexact_constructions_kernel;
using vertex_base    = CGAL::Triangulation_vertex_base_3<kernel>;
using cell_base      = CGAL::Delaunay_triangulation_cell_base_3<kernel>;
using data_structure = CGAL::Triangulation_data_structure_3<vertex_base, cell_base>;
using triangulation  = CGAL::Delaunay_triangulation_3<kernel, data_structure, CGAL::Fast_location>;
using point_type     = typename triangulation::Point;
using cell_handle    = typename triangulation::Cell_handle;
using vertex_handle  = typename triangulation::Vertex_handle;

using array_type     = Array<double, MAXVAR>;
using coord_map      = std::unordered_map<point_type, array_type>;
using value_accessor = CGAL::Data_access<coord_map>;

using LinearAlgebra = CGAL::Linear_algebraCd<double>;
using Matrix        = typename LinearAlgebra::Matrix;
using Vector        = typename LinearAlgebra::Vector;

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

    // Read
    std::stringstream file_stream;
    file_stream << file.rdbuf();
    file.close();

    // Parse Header
    std::size_t num_dim;
    std::size_t num_points;
    std::size_t num_variables;
    file_stream >> num_dim;
    if (3 != num_dim) {
        std::cerr << "ERROR: Wrong Number of Dimensions In File" << std::endl;
        std::cerr << "Expected 3 but got " << num_dim << std::endl;
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
    double x, y, z;
    xyz.resize(num_points);
    var.resize(num_points);
    for (std::size_t i = 0; i < num_points; ++i) {
        file_stream >> x;
        file_stream >> y;
        file_stream >> z;
        xyz[i] = point_type(x, y, z);
        for (std::size_t j = 0; j < num_variables; ++j) {
            file_stream >> var[i][j];
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

    // Read
    std::stringstream file_stream;
    file_stream << file.rdbuf();
    file.close();

    // Parse Header
    std::size_t num_dim;
    std::size_t num_points;
    file_stream >> num_dim;
    if (3 != num_dim) {
        std::cerr << "ERROR: Wrong Number of Dimensions In File" << std::endl;
        std::cerr << "Expected 3 but got " << num_dim << std::endl;
        std::exit(EXIT_FAILURE);
    }
    file_stream >> num_points;

    // Parse Data
    double x, y, z;
    xyz.resize(num_points);
    for (std::size_t i = 0; i < num_points; ++i) {
        file_stream >> x;
        file_stream >> y;
        file_stream >> z;
        xyz[i] = point_type(x, y, z);
    }
}

void dump_line(const point_type& xyz, const std::size_t nvar, const array_type& var) {
    for (std::size_t n = 0; n < 3; ++n) {
        std::cout << std::setw(15) << std::setprecision(8) << std::scientific << xyz[n] << " ";
    }
    for (std::size_t n = 0; n < nvar; ++n) {
        std::cout << std::setw(15) << std::setprecision(8) << std::scientific << var[n] << " ";
    }
    std::cout << std::endl;
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
        dump_line(xyz[i], nvar, var[i]);
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
    source_coords.clear();
    source_variables.clear();

    // Search for Targets
    cell_handle c;
    std::vector<array_type> target_variables(target_coords.size());
    for (std::size_t i = 0; i < target_coords.size(); ++i) {
        // Assuming some structure to the target points the
        // previous cell might be a good first guess ???
        c = dt.locate(target_coords[i], c);

        // Perform Linear Regressions
        // (X^T * X) C =  X^T y (C = Coefficients)
        //
        auto X = Matrix(4, 4);
        auto Y = Matrix(4, nvar);
        for (std::size_t j = 0; j < 4; ++j) {
            const auto p = c->vertex(j)->point();  // Get Point
            const auto v = point_to_value[p];      // Get Variables at Point

            // Copy point coords into Matrix A
            X(j, 0) = 1;
            std::copy(p.cartesian_begin(), p.cartesian_end(), X.row_begin(j) + 1);

            // Copy variables into Matrix Y
            std::copy(v.cbegin(), v.cbegin() + nvar, Y.row_begin(j));
        }
        double determinant;
        const auto Xt    = LinearAlgebra::transpose(X);
        const auto XtX   = Xt * X;
        const auto rXtX  = LinearAlgebra::inverse(XtX, determinant);
        const auto coefs = rXtX * Xt * Y;

        for (std::size_t n = 0; n < nvar; ++n) {
            target_variables[i][n] = coefs[0][n];
            for (std::size_t d = 0; d < 3; ++d) {
                target_variables[i][n] += target_coords[i][d] * coefs[d + 1][n];
            }
        }
    }

    // // Write Solution to Screen
    write_target_file_3(target_coords, nvar, target_variables);

    return EXIT_SUCCESS;
}