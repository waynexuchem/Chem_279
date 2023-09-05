    /**
     * @brief Calculate the interaction energy of a particle with its 
     * environment (all other particles in the system). 
     * 
     * @param coordinates: The coordinates of all the particles in the system as a nested list. 
     * @param i_particle: The particle index for which to calculate the energy. 
     * @param box_length: The length of the simulation box. Assume cubic box. 
     * @param cutoff: The simulation cutoff. Beyond this distance, the interaction energies are not calculated. 
     * @return double: The pariwise interaction energy of the ith particle with all other particles in the system. 
     */


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
#include "mcsim.hpp"

namespace mcsim {
    
    std::default_random_engine re;
    typedef std::vector<AtomCoord> Coordinates;
    
    // Make some types more convenient
    typedef std::array<double, 3> AtomCoord;

    /*! Generate a random double within a given range */
    double random_double(double lower_bound, double upper_bound)
    {
        std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
        return dist(re);
    }
    
    bool accept_or_reject(double delta_u, double beta){
        bool accept;
        if (delta_u <= 0) {
            accept = true;
        } else {

            // Probability that our move will be accpeted. 
            double p_acc = exp(-beta * delta_u);

            // Generate a random number on (0, 1)
            double random_number = random_double(0,1);

            if (random_number < p_acc) {
                accept = true;
            } else {
                accept = false;
            }

        }
        return accept;
    }


    /*! Generate a random integer within a given range
        The generated integer will be on the range [a,b)
    */
    int random_integer(int lower_bound, int upper_bound)
    {           
        //dist will return [a,b] but we want [a,b)
        std::uniform_int_distribution<int> dist(lower_bound, upper_bound-1);
        return dist(re);
    } 

    double pairwise_energy(double r_ij)
    {
        double r_6 = pow((1 / r_ij), 6);
        double r_12 = pow((1 / r_ij), 12);
        double pairwise_energy = 4 * (r_12 - r_6);

        return pairwise_energy;
    }

    double calculate_tail_correction(double box_length, int n_particles, double cut_off)
    {
        int N = n_particles;
        double V = box_length * box_length * box_length; 
        double rc = 1 / cut_off;
        double rc_9 = pow(rc, 9);
        double rc_3 = pow(rc, 3);
    
        double tail_correction = (8 * M_PI * N * N / (3 * V)) * ((1/3) * rc_9 - rc_3);
        
        return tail_correction;
    }


    double calculate_distance(AtomCoord coord1, AtomCoord coord2, double box_length = -1.0){
        double distance = 0;

        for (int i = 0; i < 3; i++){
            double dim_dist = abs(coord1[i] - coord2[i]);
            if (box_length > 0.0) {
                if (dim_dist > box_length) {
                    dim_dist  = remainder(dim_dist, box_length);
                }
                if (dim_dist > box_length /2) {
                    dim_dist = box_length - dim_dist;

                }
                distance += dim_dist * dim_dist;
            } else {
                distance += dim_dist * dim_dist;
            }
        }
        distance = sqrt(distance);
        return distance;
    }


    std::pair<Coordinates, double> read_xyz(std::string file_path)
    {
        // Opens up a file stream for input
        std::ifstream infile(file_path);

        // Check that it was successfully opened
        if(!infile.is_open())
        {   
            throw std::runtime_error("File path in read_xyz does not exist!");
        }
        
        double dummy; // Data that is thrown away (box length, atom indices)
        double box_length;
        int num_atoms;
        
        // Grab box_length from first number, throw the rest away
        infile >> box_length >> dummy >> dummy;
        
        // now the number of atoms
        infile >> num_atoms;
        
        // Uncomment to help troubleshoot
        //std::cout << "Box length: " << box_length << " natoms: " << num_atoms << std::endl;
        
        // Stores the atomic coordinates
        // Remember, this is a vector of arrays
        Coordinates coords;
        
        for(int i = 0; i < num_atoms; i++)
        {   
            AtomCoord coord;
            
            // Throws away the atom index
            infile >> dummy >> coord[0] >> coord[1] >> coord[2];
            
            // Add to the vector
            coords.push_back(coord);
        }

        // Makes an appropriate pair object
        return std::make_pair(coords, box_length);
    }


    double calculate_pair_energy(Coordinates coordinates, int i_particle, double box_length, double cutoff) {
    
        // int num_atoms = coordinates.size();
        AtomCoord i_position = coordinates.at(i_particle);
        double e_total = 0.0; 

        for (size_t j_particle = 0; j_particle < coordinates.size(); j_particle++) {
            if (i_particle != j_particle) {
                AtomCoord j_position = coordinates.at(j_particle);
                double rij = calculate_distance(i_position, j_position, box_length);

                if (rij < cutoff) {
                    double e_pair = pairwise_energy(rij);
                    e_total += e_pair;
                }
            }
        }

        return e_total; 
    }

    double calculate_total_energy(Coordinates coordinates, double box_length, double cutoff) {

        double total_energy = 0;
        size_t num_atoms = coordinates.size();

        for (int i = 0; i < num_atoms; i++) {
            for (int j = i + 1; j < num_atoms; j++) {
                double dist_ij = mcsim::calculate_distance(coordinates[i], coordinates[j], box_length);
                if (dist_ij < cutoff) {
                    total_energy += mcsim::pairwise_energy(dist_ij);
                }
            }
        }
        return total_energy;
    }
}
// close namespace mcsim 



