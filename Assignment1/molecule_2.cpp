// import packages
#include "molecule.h"

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

// Define a class called molecule
molecule::molecule(const std::string& file_path) {
    // Initialize num_atoms and coordinates by reading data from the XYZ file
    std::pair<Coordinates, double> result = read_xyz(file_path);
    num_atoms = result.second;
    coordinates = result.first;

} // close class molecule

// Define read_xyz function to extract coordinates from xyz file
std::pair<Coordinates, double> molecule::read_xyz(const std::string& file_path){

        // Opens up a file stream for input
        std::ifstream infile(file_path);

        // Check that it was successfully opened
        if(!infile.is_open())
        {   
            throw std::runtime_error("File path in read_xyz does not exist!");
        }
        
        int num_atoms=0;
        std::string line;
                
        // Stores the atomic coordinates as a vector of arrays
        Coordinates coords;
        AtomCoord coord;

        // Read through the file to get the first line as the number of atoms
        infile >> num_atoms;
        std::cout << "Number of atoms read: " << num_atoms << std::endl; // Print number of atoms
        int iteration = num_atoms;

        infile.clear();             // Clear the EOF flag

        // Use the acquired number of atoms to read through all the coordinates and store them in a vector
        for(int i = 0; i < iteration; i++){   

            int atom_number;

            // Throws away the atom index
            infile >> atom_number >> coord[0] >> coord[1] >> coord[2];

            // Check if the atom number equals to 79 (Au)
            if(atom_number == 79){
                // Add to the vector
                coords.push_back(coord);
            }
            else{
                std::cerr << "Non-Au atom found in the datafile!" << std::endl;
                num_atoms -= 1;
            }            
        }

        // Makes an appropriate pair object
        return std::make_pair(coords, num_atoms);
    }

// Define write_truncation_error function to write truncation errors to a file
void molecule::write_truncation_error(const std::string& file_path, std::vector<double> h_values, 
                            std::vector<double> forward_errors, std::vector<double> central_errors){

        // Opens up a file stream for output
        std::ofstream outfile(file_path);

        outfile.clear();             // Clear the EOF flag

        // Check that it was successfully opened
        if(!outfile.is_open())
        {   
            throw std::runtime_error("File path in write_truncation_error does not exist!");
        }

        // Write the data
            // Store the errors in a data file for later plotting
            // std::ofstream dataFile("truncation_errors.txt");
            // for (size_t i = 0; i < h_values.size(); i++) {
            //     dataFile << h_values[i] << " " << forward_errors[i] << " " << central_errors[i] << "\n";
            // }
            // dataFile.close();
        for (int i = 0; i < h_values.size(); i++){
            outfile << h_values[i] << "\t" << forward_errors[i] << "\t" << central_errors[i] << std::endl;
        }

    }

// Define a helper function to calculate the distance between two atoms
double molecule::calculate_distance(AtomCoord coord1, AtomCoord coord2){
        
        double r_ij = 0;

        for (int i = 0; i < 3; i++){
            double dim_dist = abs(coord1[i] - coord2[i]);
            r_ij += dim_dist * dim_dist;
            }
        r_ij = sqrt(r_ij);
        return r_ij;
    }

// Define a helper function to calculate the pairwise energy between two atoms
// given the distance between them
// Important parameters: εAu =5.29 kcal mol−1, σAu =2.951 Å
// Eij = εij * [(σij/rij)^12 − 2* (σij/rij)^6]
double molecule::pairwise_energy(double r_ij){
    double epsilon_Au = 5.29;
    double sigma_Au = 2.951;
    
    double epsilon_AuAu = sqrt(epsilon_Au * epsilon_Au);
    double sigma_AuAu = sqrt(sigma_Au * sigma_Au);

    double r_6 = pow((sigma_AuAu / r_ij), 6);
    double r_12 = pow((sigma_AuAu / r_ij), 12);
    double pairwise_energy = epsilon_AuAu * (r_12 - 2 * r_6);

    return pairwise_energy;
    }

// Define a function to calculate the total pairwise energy between all pairs of atoms
double molecule::calculate_total_energy(Coordinates coordinates, int num_atoms){

        double total_energy = 0;

        for (int i = 0; i < num_atoms; i++) {
            for (int j = i + 1; j < num_atoms; j++) {
                double r_ij = calculate_distance(coordinates[i], coordinates[j]);
                total_energy += pairwise_energy(r_ij); 
            }
        }
        return total_energy;
}

// Define a helper function to calculate the analytical force between two atoms
// given the distance between them
// Important parameters: εAu =5.29 kcal mol−1, σAu =2.951 Å
// Eij = εij * [(σij/rij)^12 − (σij/rij)^6]
// dEij/drij = 12εij * [−1/σij * (σij/rij)^13 + 1/σij * (σij/rij)^7]
double molecule::analytical_force(double r_ij){
        double epsilon_Au = 5.29;
        double sigma_Au = 2.951;
        
        double epsilon_AuAu = sqrt(epsilon_Au * epsilon_Au);
        double sigma_AuAu = sqrt(sigma_Au * sigma_Au);

        double r_7 = pow((sigma_AuAu / r_ij), 7);
        double r_13 = pow((sigma_AuAu / r_ij), 13);
        double analytical_force = 12 * epsilon_AuAu * (1/sigma_AuAu * r_7 - 1/sigma_AuAu * r_13);

        return analytical_force;
    }

// Define a function to calculate the total analytical force between all pairs of atoms
double molecule::analytical_force_total(Coordinates coordinates, int num_atoms){

        double analytical_force_total = 0;

        for (int i = 0; i < num_atoms; i++) {
            for (int j = i + 1; j < num_atoms; j++) {
                double r_ij = calculate_distance(coordinates[i], coordinates[j]);
                    analytical_force_total += analytical_force(r_ij); 
            }
        }
        return analytical_force_total;
}

// Define a helper function to calculate the forward difference force between two atoms
// given the distance between them
// Forward difference: Fi ≊h−1 [E(Ri +h) −E(Ri )]
double molecule::forward_difference(double r_ij, double h){
        
        double E_r = pairwise_energy(r_ij);
        double E_rh = pairwise_energy(r_ij + h);

        double forward_difference = 1/h * (E_rh - E_r);

        return forward_difference;
    }


// Define a function to calculate the total forward difference force between all pairs of atoms
double molecule::forward_difference_total(Coordinates coordinates, int num_atoms, double h){

        double forward_difference_total = 0;

        for (int i = 0; i < num_atoms; i++) {
            for (int j = i + 1; j < num_atoms; j++) {
                double r_ij = calculate_distance(coordinates[i], coordinates[j]);
                    forward_difference_total += forward_difference(r_ij, h);
            }
        }
        return forward_difference_total;
}

// Define a function to calculate the forward difference truncation error between all pairs of atoms
// d2Eij/drij2 = 12εij/(σij^2) * [13 * (σij/rij)^14 - 7 * (σij/rij)^8]
// et ~ 1/2 * h * d2Eij/drij2
double molecule::forward_difference_truncation_error(double r_ij, double h){
            
            double epsilon_Au = 5.29;
            double sigma_Au = 2.951;
            
            double epsilon_AuAu = sqrt(epsilon_Au * epsilon_Au);
            double sigma_AuAu = sqrt(sigma_Au * sigma_Au);
    
            double r_8 = pow((sigma_AuAu / r_ij), 8);
            double r_14 = pow((sigma_AuAu / r_ij), 14);
            double derivative_2 = 12.0 * epsilon_AuAu/(sigma_AuAu * sigma_AuAu) * (13.0 * r_14 - 7.0 * r_8);
            double forward_difference_truncation_error = 0.5 * h * derivative_2;
    
            return forward_difference_truncation_error;
}

double molecule::forward_difference_truncation_error_total(Coordinates coordinates, int num_atoms, double h){

        double forward_difference_truncation_error_total = 0;

        for (int i = 0; i < num_atoms; i++) {
            for (int j = i + 1; j < num_atoms; j++) {
                double r_ij = calculate_distance(coordinates[i], coordinates[j]);
                    forward_difference_truncation_error_total += forward_difference_truncation_error(r_ij, h);
            }
        }
        return forward_difference_truncation_error_total;
}

// Define a function to calculate the central difference force between all pairs of atoms
// Central difference: Fi ≊(2h)−1 [E(Ri +h) −E(Ri −h)]
double molecule::central_difference(double r_ij, double h){
        
        double E_rh = pairwise_energy(r_ij + h);
        double E_rm = pairwise_energy(r_ij - h);

        double central_difference = 1.0/(2.0*h) * (E_rh - E_rm);

        return central_difference;
    }

// Define a function to calculate the total central difference force between all pairs of atoms
double molecule::central_difference_total(Coordinates coordinates, int num_atoms, double h){

        double central_difference_total = 0;

        for (int i = 0; i < num_atoms; i++) {
            for (int j = i + 1; j < num_atoms; j++) {
                double r_ij = calculate_distance(coordinates[i], coordinates[j]);
                    central_difference_total += central_difference(r_ij, h);
            }
        }
        return central_difference_total;
}

// Define a function to calculate the central difference truncation error
// d3Eij/drij3 = 84εij/(σij^3) * [-26 * (σij/rij)^15 + 8 * (σij/rij)^9]
// et ~ -1/6 * h^2 * d3Eij/drij3
double molecule::central_difference_truncation_error(double r_ij, double h){
            
            double epsilon_Au = 5.29;
            double sigma_Au = 2.951;
            
            double epsilon_AuAu = sqrt(epsilon_Au * epsilon_Au);
            double sigma_AuAu = sqrt(sigma_Au * sigma_Au);
    
            double r_9 = pow((sigma_AuAu / r_ij), 9);
            double r_15 = pow((sigma_AuAu / r_ij), 15);
            double derivative_3 = 84.0 * epsilon_AuAu/(sigma_AuAu * sigma_AuAu * sigma_AuAu) * (-26.0 * r_15 + 8.0 * r_9);
            double central_difference_truncation_error = 1.0/6.0 * h * h * derivative_3;
    
            return central_difference_truncation_error;
}

double molecule::central_difference_truncation_error_total(Coordinates coordinates, int num_atoms, double h){

        double central_difference_truncation_error_total = 0;

        for (int i = 0; i < num_atoms; i++) {
            for (int j = i + 1; j < num_atoms; j++) {
                double r_ij = calculate_distance(coordinates[i], coordinates[j]);
                    central_difference_truncation_error_total += central_difference_truncation_error(r_ij, h);
            }
        }
        return central_difference_truncation_error_total;
}
    
// Define a helper function to calculate the analytical force of one individual atom
// given the coordinates of all atoms
// dE/dx_i = d/dx_i(Eij) = dEij/drij * drij/dx_i = analytical_force(double r_ij) * (x_i - x_j)/r_ij
// dE/dy_i = d/dy_i(Eij) = dEij/drij * drij/dy_i = analytical_force(double r_ij) * (y_i - y_j)/r_ij
// dE/dz_i = d/dz_i(Eij) = dEij/drij * drij/dz_i = analytical_force(double r_ij) * (z_i - z_j)/r_ij
AtomCoord molecule::analytical_force_array(AtomCoord atom, Coordinates coordinates, int num_atoms){
        AtomCoord forces = {0.0, 0.0, 0.0}; // Initialize forces to 0 for x, y, and z

        for (int i = 0; i < num_atoms; i++) {
            if (atom != coordinates[i]) { // Exclude the current atom
                double r_ij = calculate_distance(atom, coordinates[i]);
                double dx = coordinates[i][0] - atom[0];
                double dy = coordinates[i][1] - atom[1];
                double dz = coordinates[i][2] - atom[2];

                forces[0] += analytical_force(r_ij) * dx/r_ij;
                forces[1] += analytical_force(r_ij) * dy/r_ij; 
                forces[2] += analytical_force(r_ij) * dz/r_ij;
                }
            }

        // // print out the forces for debugging
        // std::cout << "Gradient forces are: ";
        // for (double value : forces) {
        //     std::cout << std::fixed << std::setprecision(6) << value << " " << std::endl;
        // }

        return forces;
        }

// Define a function to calculate the analytical force of all atoms and return a matrix
Coordinates molecule::analytical_force_matrix(Coordinates coordinates, int num_atoms){
        Coordinates forces_matrix(num_atoms, AtomCoord{0.0, 0.0, 0.0}); // Initialize forces to 0 for x, y, and z

        for (int i = 0; i < num_atoms; i++) {
            AtomCoord forces = analytical_force_array(coordinates[i], coordinates, num_atoms);
            forces_matrix[i] = forces;
        }

        // print out the forces_matrix for debugging
        std::cout << "Forces matrix is: ";
        for (const auto& coord : forces_matrix) {
            std::cout << "Matrix forces are: ";
            for (double value : coord) {
                std::cout << std::fixed << std::setprecision(6) << value << " ";
            }
            std::cout << std::endl;
        }

        return forces_matrix;
}

// Define steepest descent function to optimize gold cluster geometry and energy
Coordinates molecule::steepest_descent(Coordinates coordinates, int num_atoms, double step_size, double tolerance, int max_iterations){
        Coordinates new_coordinates = coordinates;  // Initialize new coordinates to be the same as the old ones
        
        double total_energy_old = calculate_total_energy(coordinates, num_atoms);  // Calculate total energy

        int iteration = 0;  // Initialize iteration to be 0

        while (iteration < max_iterations){

            Coordinates forces_matrix = analytical_force_matrix(coordinates, num_atoms);  // Calculate forces matrix as gradient

            for (int i = 0; i < num_atoms; i++) {
                for (int j = 0; j < 3; j++) {
                    new_coordinates[i][j] = coordinates[i][j] - step_size * forces_matrix[i][j];
                }
            }

            double total_energy_new = calculate_total_energy(new_coordinates, num_atoms);

            double total_energy_difference = std::abs(total_energy_new - total_energy_old);

            std::cout << "Iteration " << iteration << " Total energy: " << total_energy_new << std::endl;
            std::cout << "Step size: " << step_size << std::endl;
            std::cout << "Total energy difference: " << total_energy_difference << std::endl;

            if (total_energy_difference < tolerance) {
                break;  // Convergence criterion met
            }

            if (total_energy_new < total_energy_old) {
                coordinates = new_coordinates;
                step_size *= 1.1;
            } else {
                step_size /= 2.0;
                new_coordinates = coordinates;  // Revert to previous coordinates
            }

            total_energy_old = total_energy_new;
            iteration++;

        }

        return new_coordinates;
    }
