/**
 * @file run_simulation.cpp
 * @author Wang Xu, Charis Liao 
 * @brief Runs Monte Carlo simulation by num_steps times and returns the Lennard Jones Potential 
 * for each confirguration and steps as a pair. 
 * @version 0.1
 * @date 2022-08-17
 * 
 * @copyright Copyright (c) 2022
 *
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
using namespace mcsim;
using namespace std;

// Make some types more convenient
typedef array<double, 3> AtomCoord;
typedef vector<AtomCoord> Coordinates;


/**
 * @brief Runs Monte Carlo simulation to calculate the Lennard Jones energy between paritcles. 
 * 
 * @param coords: A vector of arrays of particle coordinates 
 * @param box_length: The length of the simulation box. Assume cubic box.
 * @param reduced_temperature: The simulation temperature. 
 * @param num_steps: Number of simulations.  
 * @param max_displacement: Maximum displacement for a particular particle. 
 * @param cutoff: A simulation cutoff. Beyond this distance, the interaction energies are not calculated.
 * @return vector<pair<int, double>> : A vector of <steps, particle energy> pair. 
 */

vector<pair<int, double>> run_simulation(Coordinates coords, double box_length, 
                                    double reduced_temperature, int num_steps, double max_displacement,
                                    double cutoff) {

    // Install initial parameters
    double beta = 1.0 / reduced_temperature;
    size_t num_particles = coords.size();

    // Create Empty vector to store pairs of (step, array)
    vector<pair<int, double>> result; 

    // Calculate starting energies
    double total_energy = calculate_total_energy(coords, box_length, cutoff);
    cout << "Initial total energy: " << total_energy << endl;
    double tail_correction = calculate_tail_correction(box_length, num_particles, cutoff);
    cout << "System tail correction: " << tail_correction << endl;

    // Count tail correction into total energy
    total_energy += tail_correction;

    // Loop to create system configurations
    for(size_t i = 0; i < num_steps; i++) {
                
        // 1. Randomly pick one of the particles with uniform probability
        int random_particle = random_integer(0, num_particles);

        // 2. Calculate the energy of the system
        double current_energy = calculate_pair_energy(coords, random_particle, box_length, cutoff);

        // 3. Generate a random x, y, z displacement within max_displacement
        double x_rand = random_double(-max_displacement, max_displacement);
        double y_rand = random_double(-max_displacement, max_displacement);
        double z_rand = random_double(-max_displacement, max_displacement);

        // 4. Modify the coordinates of selected particle using generate displacement
        coords.at(random_particle)[0] += x_rand;
        coords.at(random_particle)[1] += y_rand;
        coords.at(random_particle)[2] += z_rand;

        // 5. Calculate the new energy of the system
        double proposed_energy = calculate_pair_energy(coords, random_particle, box_length, cutoff);
        double delta_energy = proposed_energy - current_energy;

        // 6. Calculate if we accept the move based on Metropolis criterion
        bool accept = accept_or_reject(delta_energy, beta);

        // 7. if accepted, move the particles and update the energy
        if (accept){
            total_energy += delta_energy;
        }     
        else{
            // Move is not accepted, roll back coordinates
            coords[random_particle][0] -= x_rand;
            coords[random_particle][1] -= y_rand;
            coords[random_particle][2] -= z_rand;
        }
        if ((i+1) % 1000 == 0) {
            pair<int, double> sub_result = pair<int, double>((i+1), total_energy / num_particles);
            result.push_back(sub_result);
        }
    }
    return result;
}

int main(void) {
    std::pair<Coordinates, double> xyz_info = 
    read_xyz("lj_sample_configurations/lj_sample_config_periodic1.txt");
    Coordinates coords = xyz_info.first;
    double box_length = xyz_info.second;
    vector<pair<int, double>> sample_result = run_simulation(coords, box_length, 1, 100000, 0.1, 3);
    for (size_t i = 0; i < sample_result.size(); i++) {
        auto [step, energy] = sample_result.at(i);
        cout << "In step " << step << " the Lennard Jones Potential is: " << energy << endl;
    } 
    return 0;
}