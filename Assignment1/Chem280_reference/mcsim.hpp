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
#include <cmath>
#define _USER_MATH_DEFINES



namespace mcsim {
    
    // Make some types more convenient
    typedef std::array<double, 3> AtomCoord;
    typedef std::vector<AtomCoord> Coordinates;

    /*! Generate a random double within a given range */
    double random_double(double lower_bound, double upper_bound);
    
    bool accept_or_reject(double delta_u, double beta);

    double calculate_tail_correction(double box_length, int n_particles, double cut_off);

    double calculate_distance(AtomCoord coord1, AtomCoord coord2, double box_length);

    std::pair<Coordinates, double> read_xyz(std::string file_path);

    double pairwise_energy(double r_ij);

    double calculate_pair_energy(Coordinates coordinates, int i_particle, double box_length, double cutoff);

    double calculate_total_energy(Coordinates coordinates, double box_length, double cutoff);

    int random_integer(int lower_bound, int upper_bound);

    

}// close namespace mcsim 