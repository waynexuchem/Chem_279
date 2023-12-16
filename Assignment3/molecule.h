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

// Define quantum numbers (l, m, n) as a vector of integers
typedef std::vector<int> QuantumNumbers;

extern std::vector<double> h_exponents;
extern std::vector<double> h_coefficients;
extern std::vector<double> c_exponents;
extern std::vector<double> c_2s_coefficients;
extern std::vector<double> c_2p_coefficients;

// create a class called molecule
class molecule{
    public:
        Coordinates H_coords;
        Coordinates C_coords;

    // Constructor
        molecule(const std::string& file_path);

    // Methods
        std::pair<Coordinates, Coordinates> read_xyz(const std::string& file_path);
        int calculate_num_basis_functions(Coordinates C_coords, Coordinates H_coords);
        int calculate_electrons(Coordinates C_coords, Coordinates H_coords);
        AtomCoord calculateCenter(AtomCoord centerA, AtomCoord centerB, double alpha, double beta);
        double gaussian_func(double x, double center, double exponent, int quantumNumber);
        tuple<double, double, double> primitive_Gaussian(QuantumNumbers quantumNumbers, AtomCoord center, AtomCoord coord, vector<double> exponents);
        tuple<double, double, double> normalization_constant_S(vector<double> exponents, AtomCoord center);
        tuple<double, double, double> normalization_constant_P(vector<double> exponents, AtomCoord center);
        double calculate_basis_functions_S(vector<double> exponents, vector<double> coefficients, AtomCoord centerA, AtomCoord centerB);
        double calculate_hamiltonian_S(vector<double> exponents, vector<double> coefficients, AtomCoord centerA, AtomCoord centerB);
        std::vector<std::vector<double>> C2H2_norm_const(vector<double> h_exponents, vector<double> c_exponents, Coordinates H_coords, Coordinates C_coords);
        arma::mat calculate_basis_function_C2H2 (vector<double> h_exponents, vector<double> c_exponents, vector<double> h_coefficients, 
                                                    vector<double> c_2s_coefficients, vector<double> c_2p_coefficients, Coordinates H_coords, Coordinates C_coords);
        std::vector<std::vector<double>> C2H4_norm_const(vector<double> h_exponents, vector<double> c_exponents, Coordinates H_coords, Coordinates C_coords);
        arma::mat calculate_basis_function_C2H4 (vector<double> h_exponents, vector<double> c_exponents, vector<double> h_coefficients, 
                                                    vector<double> c_2s_coefficients, vector<double> c_2p_coefficients, Coordinates H_coords, Coordinates C_coords);
        arma::mat calculate_hamiltonian_C2H2(arma::mat overlap_matrix_C2H2);
        arma::mat calculate_hamiltonian_C2H4(arma::mat overlap_matrix_C2H4);
    }; // close class molecule