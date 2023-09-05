// header files
# pragma once

// import packages
#include <iostream> 
#include <cmath>
#include <random> // for random numbers
#include <chrono> // for generating the seed
#include <fstream> // for reading/writing files
#include <array>   // for std::array
#include <vector>  // for std::vector
#include <utility> // for std::pair
#include <string>  // for std::string
#include <iomanip>  // for printing data in tabular format
#include <armadillo>

using namespace std;

// Define individual atom coordinates as an array
typedef std::array<double, 3> AtomCoord;

// Define atom coordinates as a vector for retrieving coordinates from xyz file
typedef std::vector<AtomCoord> Coordinates;

// create a class called molecule
class molecule{
    public:
        int num_atoms;
        Coordinates coordinates;
        vector<int> atomic_numbers;

    // Constructor
        molecule(const std::string& file_path);

    // Methods
        std::pair<Coordinates, double> read_xyz(const std::string& file_path);
        void write_truncation_error(const std::string& file_path, std::vector<double> h_values, 
                        std::vector<double> forward_errors, std::vector<double> central_errors);
        double calculate_distance(AtomCoord coord1, AtomCoord coord2);
        double pairwise_energy(double r_ij);
        double calculate_total_energy(Coordinates coordinates, int num_atoms);
        double analytical_force(double r_ij);
        double analytical_force_total(Coordinates coordinates, int num_atoms);
        double forward_difference(double r_ij, double h);
        double forward_difference_total(Coordinates coordinates, int num_atoms, double h);
        double forward_difference_truncation_error(double r_ij, double h);
        double forward_difference_truncation_error_total(Coordinates coordinates, int num_atoms, double h);
        double central_difference(double r_ij, double h);
        double central_difference_total(Coordinates coordinates, int num_atoms, double h);
        double central_difference_truncation_error(double r_ij, double h);
        double central_difference_truncation_error_total(Coordinates coordinates, int num_atoms, double h);
        AtomCoord analytical_force_array(AtomCoord atom, Coordinates coordinates, int num_atoms);
        Coordinates analytical_force_matrix(Coordinates coordinates, int num_atoms);
        Coordinates steepest_descent(Coordinates coordinates, int num_atoms, double step_size, double tolerance, int max_iterations);
    }; // close class molecule
