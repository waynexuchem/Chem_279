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


#include "molecule.h"
using namespace std;

// Main function to run the program
int main(void) {

    // Input the path to the XYZ file
    std::string file_path;
    std::cout << "Enter the path to the XYZ file: ";
    std::cin >> file_path;

    // Create an instance of the molecule class
    molecule molecule(file_path);

    // Read in the coordinates and number of atoms from the XYZ file
    std::pair<Coordinates, double> xyz_info = molecule.read_xyz(file_path);
    Coordinates coords = xyz_info.first;

    // Print out the coordinates to confirm that they were read in correctly
    for (const auto& coord : coords) {
        std::cout << "Found coordinates: ";
        for (double value : coord) {
            std::cout << std::fixed << std::setprecision(6) << value << " ";
        }
        std::cout << std::endl;
    }

    // Print out the number of atoms to confirm that read_xyz worked correctly
    int num_atoms = xyz_info.second;
    std::cout << "Found number of Au atoms: " << num_atoms << std::endl;

    /**
     * Problem 1: Calculate the Lennard Jones energy of a system with multiple particles.
     * 
     * @param coordinates: The coordinates of all the particles in the system as a nested list. 
     * @return double: The total energy calculated from pariwise interaction energy of the 
     * ith particle with all other particles in the system. 
     */
    
    double total_energy = molecule.calculate_total_energy(coords, num_atoms);

    std::cout << "Problem 1 Total energy: " << total_energy << std::endl;

    /**
     * Problem 2: Calculate the analytical force from Lennard Jones energy 
     * and compare with forward / central differences.
     * 
     * @param coordinates: The coordinates of all the particles in the system as a nested list. 
     * @return double: The analytical force of any given pair of particles in the system. 
     */
    double potential_energy_surface = molecule.analytical_force_total(coords, num_atoms);
    double potential_energy_forward = molecule.forward_difference_total(coords, num_atoms, 0.01);
    double potential_energy_central = molecule.central_difference_total(coords, num_atoms, 0.01);

    std::cout << "Problem 2 Results: " << std::endl;
    std::cout << std::left << std::setw(50) << "Method" << std::setw(30) << "Potential Energy Surface" << std::endl;
    std::cout << std::left << std::setw(50) << "Analytical" << std::setw(30) << potential_energy_surface << std::endl;
    std::cout << std::left << std::setw(50) << "Forward Difference" << std::setw(30) << potential_energy_forward << std::endl;
    std::cout << std::left << std::setw(50) << "Central Difference" << std::setw(30) << potential_energy_central << std::endl;

    // prepare data for plotting in problem 2
    std::vector<double> h_values = {0.1, 0.01, 0.001, 0.0001};
    std::vector<double> forward_errors;
    std::vector<double> central_errors;

    for (double h : h_values) {
        // Calculate truncation errors using the current h value
        double forward_error = molecule.forward_difference_truncation_error_total(coords, num_atoms, h);
        double central_error = molecule.central_difference_truncation_error_total(coords, num_atoms, h);

        forward_errors.push_back(forward_error);
        central_errors.push_back(central_error);
    }

    // Print out the truncation errors to confirm that they were calculated correctly
    std::cout << "Problem 2 Truncation errors: " << std::endl;
    std::cout << std::left << std::setw(50) << "h" << std::setw(30) << "Forward Error" << std::setw(30) << "Central Error" << std::endl;
    for (int i = 0; i < h_values.size(); i++){
        std::cout << std::left << std::setw(50) << h_values[i] << std::setw(30) << forward_errors[i] << std::setw(30) << central_errors[i] << std::endl;
    }

    // Write the truncation errors to a file
    molecule.write_truncation_error("truncation_errors.csv", h_values, forward_errors, central_errors);

    /**
     * Problem 3: Use steepest descent to optimize gold cluster geometry and energy. 
     * 
     * @param coordinates: The coordinates of all the particles in the system as a nested list. 
     * @return double: The analytical force of any given pair of particles in the system. 
     */

    Coordinates force_matrix = molecule.analytical_force_matrix(coords, num_atoms);

    // Convert force_matrix to an Armadillo matrix
    arma::mat force_arma(force_matrix.size(), 3);

    for (size_t i = 0; i < force_matrix.size(); ++i) {
        for (size_t j = 0; j < 3; ++j) {
            force_arma(i, j) = force_matrix[i][j];
        }
    }

    double step_size = arma::norm(force_arma, "fro");

    std::cout << "Step_size is: " << step_size << std::endl;

    double tolerance = 1e-5;
    int max_iter = 1000;
    Coordinates new_coords = molecule.steepest_descent(coords, num_atoms, step_size, tolerance, max_iter);
    
    double new_total_energy = molecule.calculate_total_energy(new_coords, num_atoms);

    std::cout << "Problem 3 New Total energy after steepest descent optimization is: " << new_total_energy << std::endl;

    // Print out the new coordinates
    for (const auto& coord : new_coords) {
        std::cout << "New coordinates are: ";
        for (double value : coord) {
            std::cout << std::fixed << std::setprecision(6) << value << " ";
        }
        std::cout << std::endl;
    }

    // clear up memory
    coords.clear();

    return 0;
}