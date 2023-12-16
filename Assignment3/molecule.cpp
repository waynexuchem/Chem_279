// import packages
#include "molecule.h"
#include "Shell.h" // Include the Shell class header
#include "Utils.h"  // import Utils package

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

// molecule.cpp
std::vector<double> h_exponents = {3.42525091, 0.62391373, 0.16885540};
std::vector<double> h_coefficients = {0.15432897, 0.53532814, 0.44463454};

std::vector<double> c_exponents = {2.94124940, 0.68348310, 0.22228990};
std::vector<double> c_2s_coefficients = {-0.09996723, 0.39951283, 0.70011547};
std::vector<double> c_2p_coefficients = {0.15591627, 0.60768372, 0.39195739};

// Define a class called molecule
molecule::molecule(const std::string& file_path) {
    // Initialize num_atoms and coordinates by reading data from the XYZ file
    std::pair<Coordinates, Coordinates> result = read_xyz(file_path);
    H_coords = result.first;
    C_coords = result.second;

} // close class molecule

// Define read_xyz function to extract coordinates from xyz file
// format is: Element Symbol, x, y, z
// Need to differentiate carbon and hydrogen atoms and save them in different vectors
pair<Coordinates, Coordinates> molecule::read_xyz(const std::string& file_path){

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
    Coordinates H_coords;
    Coordinates C_coords;
    AtomCoord coord;

    // Read through the file to get the first line as the number of atoms
    infile >> num_atoms;

    // Ignore the remaining values on the same line
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::cout << "Number of atoms read: " << num_atoms << std::endl; // Print number of atoms
    int iteration = num_atoms;

    infile.clear();             // Clear the EOF flag

    // Use the acquired number of atoms to read through all the coordinates and store them in a vector
    for(int i = 0; i < iteration; i++){   

        int atom_number;

        // check the atom symbol and see if it is H or C
        infile >> atom_number >> coord[0] >> coord[1] >> coord[2];

        if(atom_number == 1){
            H_coords.push_back(coord);
        }
        else if(atom_number == 6){
            C_coords.push_back(coord);
        }            
    }

    // Makes an appropriate pair object
    return std::make_pair(H_coords, C_coords);
    }

// Define calculate_basis_functions function to calculate the number of basis functions
int molecule::calculate_num_basis_functions(Coordinates C_coords, Coordinates H_coords){
    
    // get the size of the vectors
    int C = C_coords.size();
    int H = H_coords.size();
    
    int num_basis_functions = 4 * C + H;
    return num_basis_functions;
    }

// Member function to calculate the number of electrons (2n)
int molecule::calculate_electrons(Coordinates C_coords, Coordinates H_coords) {
    
    // get the size of the vectors
    int C = C_coords.size();
    int H = H_coords.size();
    
    int num_electron = 4 * C + H;
    if (num_electron % 2 != 0) {
        throw std::runtime_error("Invalid molecular formula: number of total electrons must be even.");
    }
    return num_electron / 2;
    }

// Write a function to calculate center of the product shell
// R_P = (α*R_A + β*R_B) / (α+β)
AtomCoord molecule::calculateCenter(AtomCoord centerA, AtomCoord centerB, double alpha, double beta) {
    
    AtomCoord centerP;
    for (int i = 0; i < 3; ++i) {
        centerP[i] = (alpha * centerA[i] + beta * centerB[i]) / (alpha + beta);
    }

    return centerP;
    }

// Define Gaussian_type_function(x)
double molecule::gaussian_func(double x, double center, double exponent, int quantumNumber) {
    return pow(x - center, quantumNumber) * exp(-exponent * pow(x - center, 2));
    }

// Define primitive_Gaussian function
tuple<double, double, double> molecule::primitive_Gaussian(QuantumNumbers quantumNumbers, AtomCoord center, AtomCoord coord, vector<double> exponents) {
    int l = quantumNumbers[0];
    int m = quantumNumbers[1];
    int n = quantumNumbers[2];
    double x = coord[0];
    double y = coord[1];
    double z = coord[2];
    double center_x = center[0];
    double center_y = center[1];
    double center_z = center[2];
    double alpha = exponents[0];
    double beta = exponents[1];
    double gamma = exponents[2];

    // define the primitive Gaussian function
    double omega_1 = gaussian_func(x, center_x, alpha, l) * gaussian_func(y, center_y, alpha, m) * gaussian_func(z, center_z, alpha, n);
    double omega_2 = gaussian_func(x, center_x, beta, l) * gaussian_func(y, center_y, beta, m) * gaussian_func(z, center_z, beta, n);
    double omega_3 = gaussian_func(x, center_x, gamma, l) * gaussian_func(y, center_y, gamma, m) * gaussian_func(z, center_z, gamma, n);

    return make_tuple(omega_1, omega_2, omega_3);
    }

// Define normalization constant for S-type functions
tuple<double, double, double> molecule::normalization_constant_S(vector<double> exponents, AtomCoord center) {
    double alpha = exponents[0];
    double beta = exponents[1];
    double gamma = exponents[2];

    arma::vec centerA(3); // Create a 3-element arma::vec
    centerA(0) = center[0];
    centerA(1) = center[1];
    centerA(2) = center[2];
    arma::vec centerB(3); // Create a 3-element arma::vec
    centerB(0) = center[0];
    centerB(1) = center[1];
    centerB(2) = center[2];

    // Create Shell objects for s-shells A and B
    Shell shellA_1(centerA, "s", alpha);
    Shell shellB_1(centerB, "s", alpha);
    Shell shellA_2(centerA, "s", beta);
    Shell shellB_2(centerB, "s", beta);
    Shell shellA_3(centerA, "s", gamma);
    Shell shellB_3(centerB, "s", gamma);

    // Calculate the overlap integral
    double overlap_1 = calculateOverlapIntegral_ss(shellA_1, shellB_1);
    double overlap_2 = calculateOverlapIntegral_ss(shellA_2, shellB_2);
    double overlap_3 = calculateOverlapIntegral_ss(shellA_3, shellB_3);
    double norm_1 = pow(overlap_1, -0.5);
    double norm_2 = pow(overlap_2, -0.5);
    double norm_3 = pow(overlap_3, -0.5);

    return make_tuple(norm_1, norm_2, norm_3);
    }

// Define normalization constant for P-type functions
tuple<double, double, double> molecule::normalization_constant_P(vector<double> exponents, AtomCoord center) {
    double alpha = exponents[0];
    double beta = exponents[1];
    double gamma = exponents[2];

    arma::vec centerA(3); // Create a 3-element arma::vec
    centerA(0) = center[0];
    centerA(1) = center[1];
    centerA(2) = center[2];
    arma::vec centerB(3); // Create a 3-element arma::vec
    centerB(0) = center[0];
    centerB(1) = center[1];
    centerB(2) = center[2];

    // Create Shell objects for s-shells A and B
    Shell shellA_1(centerA, "p", alpha);
    Shell shellB_1(centerB, "p", alpha);
    Shell shellA_2(centerA, "p", beta);
    Shell shellB_2(centerB, "p", beta);
    Shell shellA_3(centerA, "p", gamma);
    Shell shellB_3(centerB, "p", gamma);

    // Calculate the overlap integral as a matrix of size [3, 3]
    //       px  py  pz
    // px  [[xx, xy, xz],
    // py   [yx, yy, yz],
    // pz   [zx, zy, zz]]
    arma::mat overlap_1 = calculateOverlapIntegral_pp(shellA_1, shellB_1);
    arma::mat overlap_2 = calculateOverlapIntegral_pp(shellA_2, shellB_2);
    arma::mat overlap_3 = calculateOverlapIntegral_pp(shellA_3, shellB_3);

    // 
    double norm_1 = pow(overlap_1(0, 0), -0.5);
    double norm_2 = pow(overlap_2(0, 0), -0.5);
    double norm_3 = pow(overlap_3(0, 0), -0.5);

    return make_tuple(norm_1, norm_2, norm_3);
    }

double molecule::calculate_basis_functions_S(vector<double> exponents, vector<double> coefficients, AtomCoord centerA, AtomCoord centerB){

    arma::vec centerA_arma(3); // Create a 3-element arma::vec
    centerA_arma(0) = centerA[0];
    centerA_arma(1) = centerA[1];
    centerA_arma(2) = centerA[2];
    arma::vec centerB_arma(3); // Create a 3-element arma::vec
    centerB_arma(0) = centerB[0];
    centerB_arma(1) = centerB[1];
    centerB_arma(2) = centerB[2];

    // First get center of the two H atoms then get normalization constants
    const AtomCoord center = calculateCenter(centerA, centerB, exponents[0], exponents[0]);
    // create a vector to store the normalization constants
    const auto [norm_1, norm_2, norm_3] = normalization_constant_S(exponents, center);
    vector<double> norm = {norm_1, norm_2, norm_3};

    double norm_off_diag = 0.0;
    // Calculate the basis functions
    for (int k=0; k < 3; k++){
        for (int l=0; l<3; l++){
            Shell shellA(centerA_arma, "s", exponents[k]);
            Shell shellB(centerB_arma, "s", exponents[l]);
            norm_off_diag += coefficients[k] * coefficients[l] * norm[k] * norm[l] * calculateOverlapIntegral_ss(shellA, shellB);
        }
    }
    return norm_off_diag;
}

// calculate Hamiltonian matrix elements
double molecule::calculate_hamiltonian_S(vector<double> exponents, vector<double> coefficients, AtomCoord centerA, AtomCoord centerB) {
    double K = 1.75;
    double diag_Hamiltonian = -13.6;  // for H atom only

    double norm_off_diag = calculate_basis_functions_S(exponents, coefficients, centerA, centerB);
    double hamiltonian_off_diag = 0.5 * K * (diag_Hamiltonian + diag_Hamiltonian) * norm_off_diag;

    return hamiltonian_off_diag;
}

// Calculate normalization constants for C2H2 as a 10*3 matrix
// rows go as: [Hs, Cs, Cpx, Cpy, Cpz, Cs, Cpx, Cpy, Cpz, Hs]
// columns go as: [norm_1, norm_2, norm_3]
std::vector<std::vector<double>> molecule::C2H2_norm_const(vector<double> h_exponents, vector<double> c_exponents, Coordinates H_coords, Coordinates C_coords) {
    std::vector<std::vector<double>> C2H2_norm_const;

    // H_coords has 2 AtomCoord objects, C_coords has 2 AtomCoord objects
    // H_coords[0] is the first H atom, H_coords[1] is the second H atom
    // C_coords[0] is the first C atom, C_coords[1] is the second C atom
    // Sequence of calculation: H_coords[0], C_coords[0], C_coords[1], H_coords[1]
    // H_coords[0] yields 1 tuple of norm_consts, C_coords[0] yields 4 tuples of norm_consts, C_coords[1] yields 4 tuples of norm_consts, H_coords[1] yields 1 tuple of norm_consts

    // First get center of H_coords[0] then get normalization constants
    const AtomCoord center_H1 = calculateCenter(H_coords[0], H_coords[0], h_exponents[0], h_exponents[0]);
    const auto [norm_H1_1, norm_H1_2, norm_H1_3] = normalization_constant_S(h_exponents, center_H1);
    vector<double> norm_H1 = {norm_H1_1, norm_H1_2, norm_H1_3};
    C2H2_norm_const.push_back(norm_H1);

    // First get center of C_coords[0] then get normalization constants
    const AtomCoord center_Cs1 = calculateCenter(C_coords[0], C_coords[0], c_exponents[0], c_exponents[0]);
    const auto [norm_Cs1_1, norm_Cs1_2, norm_Cs1_3] = normalization_constant_S(c_exponents, center_Cs1);
    vector<double> norm_Cs1 = {norm_Cs1_1, norm_Cs1_2, norm_Cs1_3};
    C2H2_norm_const.push_back(norm_Cs1);

    const AtomCoord center_Cpx1 = calculateCenter(C_coords[0], C_coords[0], c_exponents[0], c_exponents[0]);
    const auto [norm_Cpx1_1, norm_Cpx1_2, norm_Cpx1_3] = normalization_constant_P(c_exponents, center_Cpx1);
    vector<double> norm_Cpx1 = {norm_Cpx1_1, norm_Cpx1_2, norm_Cpx1_3};
    C2H2_norm_const.push_back(norm_Cpx1);

    const AtomCoord center_Cpy1 = calculateCenter(C_coords[0], C_coords[0], c_exponents[1], c_exponents[1]);
    const auto [norm_Cpy1_1, norm_Cpy1_2, norm_Cpy1_3] = normalization_constant_P(c_exponents, center_Cpy1);
    vector<double> norm_Cpy1 = {norm_Cpy1_1, norm_Cpy1_2, norm_Cpy1_3};
    C2H2_norm_const.push_back(norm_Cpy1);

    const AtomCoord center_Cpz1 = calculateCenter(C_coords[0], C_coords[0], c_exponents[2], c_exponents[2]);
    const auto [norm_Cpz1_1, norm_Cpz1_2, norm_Cpz1_3] = normalization_constant_P(c_exponents, center_Cpz1);
    vector<double> norm_Cpz1 = {norm_Cpz1_1, norm_Cpz1_2, norm_Cpz1_3};
    C2H2_norm_const.push_back(norm_Cpz1);

    // First get center of C_coords[1] then get normalization constants
    const AtomCoord center_Cs2 = calculateCenter(C_coords[1], C_coords[1], c_exponents[0], c_exponents[0]);
    const auto [norm_Cs2_1, norm_Cs2_2, norm_Cs2_3] = normalization_constant_S(c_exponents, center_Cs2);
    vector<double> norm_Cs2 = {norm_Cs2_1, norm_Cs2_2, norm_Cs2_3};
    C2H2_norm_const.push_back(norm_Cs2);

    const AtomCoord center_Cpx2 = calculateCenter(C_coords[1], C_coords[1], c_exponents[0], c_exponents[0]);
    const auto [norm_Cpx2_1, norm_Cpx2_2, norm_Cpx2_3] = normalization_constant_P(c_exponents, center_Cpx2);
    vector<double> norm_Cpx2 = {norm_Cpx2_1, norm_Cpx2_2, norm_Cpx2_3};
    C2H2_norm_const.push_back(norm_Cpx2);

    const AtomCoord center_Cpy2 = calculateCenter(C_coords[1], C_coords[1], c_exponents[1], c_exponents[1]);
    const auto [norm_Cpy2_1, norm_Cpy2_2, norm_Cpy2_3] = normalization_constant_P(c_exponents, center_Cpy2);
    vector<double> norm_Cpy2 = {norm_Cpy2_1, norm_Cpy2_2, norm_Cpy2_3};
    C2H2_norm_const.push_back(norm_Cpy2);

    const AtomCoord center_Cpz2 = calculateCenter(C_coords[1], C_coords[1], c_exponents[2], c_exponents[2]);
    const auto [norm_Cpz2_1, norm_Cpz2_2, norm_Cpz2_3] = normalization_constant_P(c_exponents, center_Cpz2);
    vector<double> norm_Cpz2 = {norm_Cpz2_1, norm_Cpz2_2, norm_Cpz2_3};
    C2H2_norm_const.push_back(norm_Cpz2);

    // First get center of H_coords[1] then get normalization constants
    const AtomCoord center_H2 = calculateCenter(H_coords[1], H_coords[1], h_exponents[0], h_exponents[0]);
    const auto [norm_H2_1, norm_H2_2, norm_H2_3] = normalization_constant_S(h_exponents, center_H2);
    vector<double> norm_H2 = {norm_H2_1, norm_H2_2, norm_H2_3};
    C2H2_norm_const.push_back(norm_H2);

    return C2H2_norm_const;

}
// double calculateOverlapIntegral_1D(int dimension, arma::vec centerA, arma::vec centerB, double alpha, double beta, int lA, int lB)
// Function to calculate 3d overlap integral
double calculateOverlapIntegral_3D(arma::vec centerA, arma::vec centerB, vector<double> exponents_alpha, vector<double> exponents_beta, int lA, int lB, int k, int l) {
    
    double overlap_integral = calculateOverlapIntegral_1D(0, centerA, centerB, exponents_alpha[k], exponents_beta[l], lA, lB) \
                            * calculateOverlapIntegral_1D(1, centerA, centerB, exponents_alpha[k], exponents_beta[l], lA, lB) \
                            * calculateOverlapIntegral_1D(2, centerA, centerB, exponents_alpha[k], exponents_beta[l], lA, lB);
    return overlap_integral;
}

// calculate the basis functions of C2H2 as a 10*10 matrix
// rows go as: [Hs, Cs, Cpx, Cpy, Cpz, Cs, Cpx, Cpy, Cpz, Hs]
// columns go as: [Hs, Cs, Cpx, Cpy, Cpz, Cs, Cpx, Cpy, Cpz, Hs]
arma::mat molecule::calculate_basis_function_C2H2 (vector<double> h_exponents, vector<double> c_exponents, vector<double> h_coefficients, 
                                                    vector<double> c_2s_coefficients, vector<double> c_2p_coefficients, Coordinates H_coords, Coordinates C_coords){
    
    // create a 10*10 matrix to store the basis functions
    arma::mat basis_functions = arma::zeros<arma::mat>(10, 10);
    const std::vector<std::vector<double>> norm_const = C2H2_norm_const(h_exponents, c_exponents, H_coords, C_coords);
    arma::vec centerA_arma(3); // Create a 3-element arma::vec
    arma::vec centerB_arma(3); // Create a 3-element arma::vec
    string shell_symbol_A;
    string shell_symbol_B;
    vector<double> exponent_A;
    vector<double> exponent_B;
    vector<double> coefficient_A;
    vector<double> coefficient_B;
    int a;
    int b;

    // loop through the rows and columns to calculate the basis functions
    // starting from the first row: HsCs, HsCpx, HsCpy, HsCpz, HsCs, HsCpx, HsCpy, HsCpz, HsHs, 9 elements in total
    // then the second row: CsCpx, CsCpy, CsCpz, CsCs, CsCpx, CsCpy, CsCpz, CsHs, 8 elements in total
    // then the third row: CpxCpy, CpxCpz, CpxCs, CpxCpx, CpxCpy, CpxCpz, CpxHs, 7 elements in total, so on and so forth
    // the last row is only CpzsHs, 1 element in total, so outer loop i is 9 times, from 9 to 1
    // the inner loop j is from i+1 to 9
    // within each loop of i and j, need to differentiate between Hs, Cs, Cpx, Cpy, Cpz
    // if i=0, exponents/coefficients of atom 1 is chosen from h_exponents/h_coefficients, shell is S
    // if i=1 or i=5, exponents/coefficients of atom 1 is chosen from c_exponents/c_2s_coefficients, shell is s
    // if i=2 or i=3 or i=4 or i=6 or i=7 or i=8, exponents/coefficients of atom 1 is chosen from c_exponents/c_2p_coefficients, shell is p
    // if j=1 or j=5, exponents/coefficients of atom 2 is chosen from c_exponents/c_2s_coefficients, shell is s
    // if j=2 or j=3 or j=4 or j=6 or j=7 or j=8, exponents/coefficients of atom 2 is chosen from c_exponents/c_2p_coefficients, shell is p
    // if j=9, exponents/coefficients of atom 2 is chosen from h_exponents/h_coefficients, shell is s

    for (int i = 0; i < 9; i++){ // define parameters of the rows
        for (int j = i+1; j < 10; j++){
            if (i == 0){
                exponent_A = h_exponents;
                coefficient_A = h_coefficients;
                centerA_arma(0) = H_coords[0][0];
                centerA_arma(1) = H_coords[0][1];
                centerA_arma(2) = H_coords[0][2];
                shell_symbol_A = "s";
                a = 0;
                } else if (i == 1){
                    exponent_A = c_exponents;
                    coefficient_A = c_2s_coefficients;
                    centerA_arma(0) = C_coords[0][0];
                    centerA_arma(1) = C_coords[0][1];
                    centerA_arma(2) = C_coords[0][2];
                    shell_symbol_A = "s";
                    a = 0;
                    } else if (i == 2){
                        exponent_A = c_exponents;
                        coefficient_A = c_2p_coefficients;
                        centerA_arma(0) = C_coords[0][0];
                        centerA_arma(1) = C_coords[0][1];
                        centerA_arma(2) = C_coords[0][2];
                        shell_symbol_A = "p";
                        a = 0;  // px
                        
                        } else if (i == 3){
                        exponent_A = c_exponents;
                        coefficient_A = c_2p_coefficients;
                        centerA_arma(0) = C_coords[0][0];
                        centerA_arma(1) = C_coords[0][1];
                        centerA_arma(2) = C_coords[0][2];
                        shell_symbol_A = "p";
                        a = 1; // py
                        
                        } else if (i == 4){
                        exponent_A = c_exponents;
                        coefficient_A = c_2p_coefficients;
                        centerA_arma(0) = C_coords[0][0];
                        centerA_arma(1) = C_coords[0][1];
                        centerA_arma(2) = C_coords[0][2];
                        shell_symbol_A = "p";
                        a = 2; // pz
                        
                        } else if (i == 5) {
                            exponent_A = c_exponents;
                            coefficient_A = c_2s_coefficients;
                            centerA_arma(0) = C_coords[1][0];
                            centerA_arma(1) = C_coords[1][1];
                            centerA_arma(2) = C_coords[1][2];
                            shell_symbol_A = "s";
                            
                            } else if (i == 6) {
                                exponent_A = c_exponents;
                                coefficient_A = c_2p_coefficients;
                                centerA_arma(0) = C_coords[1][0];
                                centerA_arma(1) = C_coords[1][1];
                                centerA_arma(2) = C_coords[1][2];
                                shell_symbol_A = "p";
                                a = 0; // px
                                
                                } else if (i == 7) {
                                exponent_A = c_exponents;
                                coefficient_A = c_2p_coefficients;
                                centerA_arma(0) = C_coords[1][0];
                                centerA_arma(1) = C_coords[1][1];
                                centerA_arma(2) = C_coords[1][2];
                                shell_symbol_A = "p";
                                a = 1; // py
                                } else if (i == 8) {
                                exponent_A = c_exponents;
                                coefficient_A = c_2p_coefficients;
                                centerA_arma(0) = C_coords[1][0];
                                centerA_arma(1) = C_coords[1][1];
                                centerA_arma(2) = C_coords[1][2];
                                shell_symbol_A = "p";
                                a = 2; // pz
                                }
                            
            if (j == 1) { 
                exponent_B = c_exponents;
                coefficient_B = c_2s_coefficients;
                centerB_arma(0) = C_coords[0][0];
                centerB_arma(1) = C_coords[0][1];
                centerB_arma(2) = C_coords[0][2];
                shell_symbol_B = "s";
                } else if (j == 2) {
                    exponent_B = c_exponents;
                    coefficient_B = c_2p_coefficients;
                    centerB_arma(0) = C_coords[0][0];
                    centerB_arma(1) = C_coords[0][1];
                    centerB_arma(2) = C_coords[0][2];
                    shell_symbol_B = "p";
                    b = 0; // px
                    } else if (j == 3) {
                    exponent_B = c_exponents;
                    coefficient_B = c_2p_coefficients;
                    centerB_arma(0) = C_coords[0][0];
                    centerB_arma(1) = C_coords[0][1];
                    centerB_arma(2) = C_coords[0][2];
                    shell_symbol_B = "p";
                    b = 1; // py
                    } else if (j == 4) {
                    exponent_B = c_exponents;
                    coefficient_B = c_2p_coefficients;
                    centerB_arma(0) = C_coords[0][0];
                    centerB_arma(1) = C_coords[0][1];
                    centerB_arma(2) = C_coords[0][2];
                    shell_symbol_B = "p";
                    b = 2; // pz
                    } else if (j == 5) {
                        exponent_B = c_exponents;
                        coefficient_B = c_2s_coefficients;
                        centerB_arma(0) = C_coords[1][0];
                        centerB_arma(1) = C_coords[1][1];
                        centerB_arma(2) = C_coords[1][2];
                        shell_symbol_B = "s";
                        } else if (j == 6) {
                            exponent_B = c_exponents;
                            coefficient_B = c_2p_coefficients;
                            centerB_arma(0) = C_coords[1][0];
                            centerB_arma(1) = C_coords[1][1];
                            centerB_arma(2) = C_coords[1][2];
                            shell_symbol_B = "p";
                            b = 0; // px
                            } else if (j == 7) {
                            exponent_B = c_exponents;
                            coefficient_B = c_2p_coefficients;
                            centerB_arma(0) = C_coords[1][0];
                            centerB_arma(1) = C_coords[1][1];
                            centerB_arma(2) = C_coords[1][2];
                            shell_symbol_B = "p";
                            b = 1; // py
                            } else if (j == 8) {
                            exponent_B = c_exponents;
                            coefficient_B = c_2p_coefficients;
                            centerB_arma(0) = C_coords[1][0];
                            centerB_arma(1) = C_coords[1][1];
                            centerB_arma(2) = C_coords[1][2];
                            shell_symbol_B = "p";
                            b = 2; // pz
                            } else if (j == 9) {
                                exponent_B = h_exponents;
                                coefficient_B = h_coefficients;
                                centerB_arma(0) = H_coords[1][0];
                                centerB_arma(1) = H_coords[1][1];
                                centerB_arma(2) = H_coords[1][2];
                                shell_symbol_B = "s";
                                } 

            // within the two outer loops, there are two additional loops of k and l, each loop has 3 elements
            // general equationto calculate the basis functions: ∑∑h_coeff * c_coeff * norm_const * norm_const * overlap_integral
            // overlap_integral calculation depends on Shell definition
            double basis_function_ij = 0.0;
            // Calculate the basis functions

            // cout << "exponent_A: " << exponent_A[0] << "  " << exponent_A[1] << "  " << exponent_A[2] << endl;
            // cout << "exponent_B: " << exponent_B[0] << "  " << exponent_B[1] << "  " << exponent_B[2] << endl;
            // cout << "coefficient_A: " << coefficient_A[0] << "  " << coefficient_A[1] << "  " << coefficient_A[2] << endl;
            // cout << "coefficient_B: " << coefficient_B[0] << "  " << coefficient_B[1] << "  " << coefficient_B[2] << endl;

            for (int k=0; k < 3; k++){
                for (int l=0; l < 3; l++){

                    // cout << "k: " << k << "  l: " << l << endl;

                    Shell shellA(centerA_arma, shell_symbol_A, exponent_A[k]);
                    Shell shellB(centerB_arma, shell_symbol_B, exponent_B[l]);

                    // Get the quantum numbers of the two shells
                    arma::imat quantumNumberA = shellA.getQuantumNumber();
                    arma::imat quantumNumberB = shellB.getQuantumNumber();

                    // double overlap = calculateOverlapIntegral_3D(centerA_arma, centerB_arma, exponent_A, exponent_B, quantumNumberA(k, l), quantumNumberB(k, l), k, l);

                    
                    if (shell_symbol_A == "s" && shell_symbol_B == "s") {
                        double overlap = calculateOverlapIntegral_ss(shellA, shellB);
                        basis_function_ij += coefficient_A[k] * coefficient_B[l] * norm_const[i][k] * norm_const[i][l] * overlap;
                    } else if (shell_symbol_A == "p" && shell_symbol_B == "p") {
                        arma::mat overlap = calculateOverlapIntegral_pp(shellA, shellB);
                        basis_function_ij += coefficient_A[k] * coefficient_B[l] * norm_const[i][k] * norm_const[i][l] * overlap(a, b);
                    } else if (shell_symbol_A == "s" && shell_symbol_B == "p") {
                        std::vector<double> overlap = calculateOverlapIntegral_sp(shellA, shellB);
                        basis_function_ij += coefficient_A[k] * coefficient_B[l] * norm_const[i][k] * norm_const[i][l] * overlap[b];
                    } else if (shell_symbol_A == "p" && shell_symbol_B == "s") {
                        std::vector<double> overlap = calculateOverlapIntegral_sp(shellA, shellB);
                        basis_function_ij += coefficient_A[k] * coefficient_B[l] * norm_const[i][k] * norm_const[i][l] * overlap[a];
                    }
                }
            }
            // fill in matrix elements of basis_function_ij
            basis_functions(i, j) = basis_function_ij;
            basis_functions(j, i) = basis_function_ij;
            }  
        }

    // define when i = j, basis function is 1
    for (int i=0; i<10; i++){
        basis_functions(i, i) = 1.0;
    }

    return basis_functions;
}

// calculate the hamiltonian of C2H2 as a 10*10 matrix
// rows go as: [Hs, Cs, Cpx, Cpy, Cpz, Cs, Cpx, Cpy, Cpz, Hs]
// columns go as: [Hs, Cs, Cpx, Cpy, Cpz, Cs, Cpx, Cpy, Cpz, Hs]
    // ({"H_s", -13.6});
    // ({"C_s", -21.4});
    // ({"C_px", -11.4});
    // ({"C_py", -11.4});
    // ({"C_pz", -11.4});
    // H(i, j) = K / 2 * (H_i + H_j) * S(i, j)
arma::mat molecule::calculate_hamiltonian_C2H2(arma::mat overlap_matrix_C2H2){
    
    // create a 10*10 matrix to store the hamiltonians
    arma::mat hamiltonian = arma::zeros<arma::mat>(10, 10);
    double K = 1.75;
    double h_1;
    double h_2;

    for (int i=0; i<9; i++){ // define parameters of the rows
        for (int j = i+1; j < 10; j++){
            if (i == 0){
                h_1 = -13.6;
                    
                } else if (i == 1){
                    h_1 = -21.4;
                    
                    } else if (i == 2 || i == 3 || i == 4){
                        h_1 = -11.4;
                        
                        } else if (i == 5) {
                            h_1 = -21.4;
                            
                            } else if (i == 6 || i == 7 || i == 8) {
                                h_1 = -11.4;
                                
                                }
                            
            if (j == 1){ // define parameters of the columns
                h_2 = -21.4;
                } else if (j == 2 || j == 3 || j == 4) {
                    h_2 = -11.4;
                    } else if (j == 5) {
                        h_2 = -21.4;
                        } else if (j == 6 || j == 7 || j == 8) {
                            h_2 = -11.4;
                            } else if (j == 9) {
                                h_2 = -13.6;
                                }

            // fill in matrix elements of basis_function_ij
            hamiltonian(i, j) = K / 2 * (h_1 + h_2) * overlap_matrix_C2H2(i, j);
            hamiltonian(j, i) = K / 2 * (h_1 + h_2) * overlap_matrix_C2H2(i, j);
            }
        }

        // define when i = j, basis function is 1
        for (int i=0; i<10; i++){
            hamiltonian(i, i) = h_1;
        }

        return hamiltonian;
}

// Calculate normalization constants for C2H4 as a 12*3 matrix
// rows go as: [Hs, Hs, Cs, Cpx, Cpy, Cpz, Cs, Cpx, Cpy, Cpz, Hs, Hs]
// columns go as: [norm_1, norm_2, norm_3]
std::vector<std::vector<double>> molecule::C2H4_norm_const(vector<double> h_exponents, vector<double> c_exponents, Coordinates H_coords, Coordinates C_coords) {
    std::vector<std::vector<double>> C2H4_norm_const;

    // H_coords has 4 AtomCoord objects, C_coords has 2 AtomCoord objects
    // H_coords[0] is the first H atom, H_coords[1] is the second H atom
    // H_coords[2] is the third H atom, H_coords[3] is the fourth H atom
    // C_coords[0] is the first C atom, C_coords[1] is the second C atom
    // Sequence of calculation: H_coords[0], H_coords[1], C_coords[0], C_coords[1], H_coords[2], H_coords[3]
    // H_coords[0] yields 1 tuple of norm_consts, H_coords[1] yields 1 tuple of norm_consts, C_coords[0] yields 4 tuples of norm_consts, 
    // C_coords[1] yields 4 tuples of norm_consts, H_coords[2] yields 1 tuple of norm_consts, H_coords[3] yields 1 tuple of norm_consts, 

    // First get center of H_coords[0] then get normalization constants
    const AtomCoord center_H1 = calculateCenter(H_coords[0], H_coords[0], h_exponents[0], h_exponents[0]);
    const auto [norm_H1_1, norm_H1_2, norm_H1_3] = normalization_constant_S(h_exponents, center_H1);
    vector<double> norm_H1 = {norm_H1_1, norm_H1_2, norm_H1_3};
    C2H4_norm_const.push_back(norm_H1);

    const AtomCoord center_H2 = calculateCenter(H_coords[1], H_coords[1], h_exponents[0], h_exponents[0]);
    const auto [norm_H2_1, norm_H2_2, norm_H2_3] = normalization_constant_S(h_exponents, center_H2);
    vector<double> norm_H2 = {norm_H2_1, norm_H2_2, norm_H2_3};
    C2H4_norm_const.push_back(norm_H2);

    // First get center of C_coords[0] then get normalization constants
    const AtomCoord center_Cs1 = calculateCenter(C_coords[0], C_coords[0], c_exponents[0], c_exponents[0]);
    const auto [norm_Cs1_1, norm_Cs1_2, norm_Cs1_3] = normalization_constant_S(c_exponents, center_Cs1);
    vector<double> norm_Cs1 = {norm_Cs1_1, norm_Cs1_2, norm_Cs1_3};
    C2H4_norm_const.push_back(norm_Cs1);

    const AtomCoord center_Cpx1 = calculateCenter(C_coords[0], C_coords[0], c_exponents[0], c_exponents[0]);
    const auto [norm_Cpx1_1, norm_Cpx1_2, norm_Cpx1_3] = normalization_constant_P(c_exponents, center_Cpx1);
    vector<double> norm_Cpx1 = {norm_Cpx1_1, norm_Cpx1_2, norm_Cpx1_3};
    C2H4_norm_const.push_back(norm_Cpx1);

    const AtomCoord center_Cpy1 = calculateCenter(C_coords[0], C_coords[0], c_exponents[1], c_exponents[1]);
    const auto [norm_Cpy1_1, norm_Cpy1_2, norm_Cpy1_3] = normalization_constant_P(c_exponents, center_Cpy1);
    vector<double> norm_Cpy1 = {norm_Cpy1_1, norm_Cpy1_2, norm_Cpy1_3};
    C2H4_norm_const.push_back(norm_Cpy1);

    const AtomCoord center_Cpz1 = calculateCenter(C_coords[0], C_coords[0], c_exponents[2], c_exponents[2]);
    const auto [norm_Cpz1_1, norm_Cpz1_2, norm_Cpz1_3] = normalization_constant_P(c_exponents, center_Cpz1);
    vector<double> norm_Cpz1 = {norm_Cpz1_1, norm_Cpz1_2, norm_Cpz1_3};
    C2H4_norm_const.push_back(norm_Cpz1);

    // First get center of C_coords[1] then get normalization constants
    const AtomCoord center_Cs2 = calculateCenter(C_coords[1], C_coords[1], c_exponents[0], c_exponents[0]);
    const auto [norm_Cs2_1, norm_Cs2_2, norm_Cs2_3] = normalization_constant_S(c_exponents, center_Cs2);
    vector<double> norm_Cs2 = {norm_Cs2_1, norm_Cs2_2, norm_Cs2_3};
    C2H4_norm_const.push_back(norm_Cs2);

    const AtomCoord center_Cpx2 = calculateCenter(C_coords[1], C_coords[1], c_exponents[0], c_exponents[0]);
    const auto [norm_Cpx2_1, norm_Cpx2_2, norm_Cpx2_3] = normalization_constant_P(c_exponents, center_Cpx2);
    vector<double> norm_Cpx2 = {norm_Cpx2_1, norm_Cpx2_2, norm_Cpx2_3};
    C2H4_norm_const.push_back(norm_Cpx2);

    const AtomCoord center_Cpy2 = calculateCenter(C_coords[1], C_coords[1], c_exponents[1], c_exponents[1]);
    const auto [norm_Cpy2_1, norm_Cpy2_2, norm_Cpy2_3] = normalization_constant_P(c_exponents, center_Cpy2);
    vector<double> norm_Cpy2 = {norm_Cpy2_1, norm_Cpy2_2, norm_Cpy2_3};
    C2H4_norm_const.push_back(norm_Cpy2);

    const AtomCoord center_Cpz2 = calculateCenter(C_coords[1], C_coords[1], c_exponents[2], c_exponents[2]);
    const auto [norm_Cpz2_1, norm_Cpz2_2, norm_Cpz2_3] = normalization_constant_P(c_exponents, center_Cpz2);
    vector<double> norm_Cpz2 = {norm_Cpz2_1, norm_Cpz2_2, norm_Cpz2_3};
    C2H4_norm_const.push_back(norm_Cpz2);

    // First get center of H_coords[2] then get normalization constants
    const AtomCoord center_H3 = calculateCenter(H_coords[2], H_coords[2], h_exponents[0], h_exponents[0]);
    const auto [norm_H3_1, norm_H3_2, norm_H3_3] = normalization_constant_S(h_exponents, center_H3);
    vector<double> norm_H3 = {norm_H3_1, norm_H3_2, norm_H3_3};
    C2H4_norm_const.push_back(norm_H3);

    // First get center of H_coords[3] then get normalization constants
    const AtomCoord center_H4 = calculateCenter(H_coords[3], H_coords[3], h_exponents[0], h_exponents[0]);
    const auto [norm_H4_1, norm_H4_2, norm_H4_3] = normalization_constant_S(h_exponents, center_H4);
    vector<double> norm_H4 = {norm_H4_1, norm_H4_2, norm_H4_3};
    C2H4_norm_const.push_back(norm_H4);

    return C2H4_norm_const;

}

// calculate the basis functions of C2H4 as a 12*12 matrix
// rows go as: [Hs, Hs, Cs, Cpx, Cpy, Cpz, Cs, Cpx, Cpy, Cpz, Hs, Hs]
// columns go as: [Hs, Hs, Cs, Cpx, Cpy, Cpz, Cs, Cpx, Cpy, Cpz, Hs, Hs]
arma::mat molecule::calculate_basis_function_C2H4 (vector<double> h_exponents, vector<double> c_exponents, vector<double> h_coefficients, 
                                                    vector<double> c_2s_coefficients, vector<double> c_2p_coefficients, Coordinates H_coords, Coordinates C_coords){
    
    // create a 12*12 matrix to store the basis functions
    arma::mat basis_functions = arma::zeros<arma::mat>(12, 12);
    const std::vector<std::vector<double>> norm_const = C2H4_norm_const(h_exponents, c_exponents, H_coords, C_coords);
    arma::vec centerA_arma(3); // Create a 3-element arma::vec
    arma::vec centerB_arma(3); // Create a 3-element arma::vec
    string shell_symbol_A;
    string shell_symbol_B;
    vector<double> exponent_A;
    vector<double> exponent_B;
    vector<double> coefficient_A;
    vector<double> coefficient_B;
    int a;
    int b;

    // loop through the rows and columns to calculate the basis functions
    // starting from the first row: HsHs, HsCs, HsCpx, HsCpy, HsCpz, HsCs, HsCpx, HsCpy, HsCpz, HsHs, HsHs, 11 elements in total
    // then the second row: HsCs, HsCpx, HsCpy, HsCpz, HsCs, HsCpx, HsCpy, HsCpz, HsHs, HsHs, 10 elements in total
    // then the third row: CsCpx, CsCpy, CsCpz, CsCs, CsCpx, CsCpy, CsCpz, CsHs, CsHs, 9 elements in total, so on and so forth
    // the last row is only HsHs, 1 element in total, so outer loop i is 11 times, from 0 to 10
    // the inner loop j is from i+1 to 11
    // within each loop of i and j, need to differentiate between Hs, Cs, Cpx, Cpy, Cpz
    // if i=0, exponents/coefficients of atom 1 is chosen from h_exponents/h_coefficients, shell is S, H_coords[0]
    // if i=1, exponents/coefficients of atom 1 is chosen from h_exponents/h_coefficients, shell is S, H_coords[1]
    // if i=2, exponents/coefficients of atom 1 is chosen from c_exponents/c_2s_coefficients, shell is s, C_coords[0]
    // if i=3 or i=4 or i=5 exponents/coefficients of atom 1 is chosen from c_exponents/c_2p_coefficients, shell is p, C_coords[0]
    // if i=6, exponents/coefficients of atom 1 is chosen from c_exponents/c_2s_coefficients, shell is s, C_coords[1]
    // if i=7 or i=8 or i=9 exponents/coefficients of atom 1 is chosen from c_exponents/c_2p_coefficients, shell is p, C_coords[1]
    // if i=10, exponents/coefficients of atom 1 is chosen from h_exponents/h_coefficients, shell is S, H_coords[2]
    // if j=1, exponents/coefficients of atom 2 is chosen from h_exponents/h_coefficients, shell is S, H_coords[1]
    // if j=2, exponents/coefficients of atom 2 is chosen from c_exponents/c_2s_coefficients, shell is s, C_coords[0]
    // if j=3 or j=4 or j=5, exponents/coefficients of atom 2 is chosen from c_exponents/c_2p_coefficients, shell is p, C_coords[0]
    // if j=6, exponents/coefficients of atom 2 is chosen from c_exponents/c_2s_coefficients, shell is s, C_coords[1]
    // if j=7 or j=8 or j=9, exponents/coefficients of atom 2 is chosen from c_exponents/c_2p_coefficients, shell is p, C_coords[1]
    // if j=10, exponents/coefficients of atom 2 is chosen from h_exponents/h_coefficients, shell is S, H_coords[2]
    // if j=11, exponents/coefficients of atom 2 is chosen from h_exponents/h_coefficients, shell is S, H_coords[3]

    for (int i=0; i<11; i++){
        for (int j = i+1; j < 12; j++){
            if (i == 0){
                exponent_A = h_exponents;
                coefficient_A = h_coefficients;
                centerA_arma(0) = H_coords[0][0];
                centerA_arma(1) = H_coords[0][1];
                centerA_arma(2) = H_coords[0][2];
                shell_symbol_A = "s";
                    
                } else if (i == 1) {
                    exponent_A = h_exponents;
                    coefficient_A = h_coefficients;
                    centerA_arma(0) = H_coords[1][0];
                    centerA_arma(1) = H_coords[1][1];
                    centerA_arma(2) = H_coords[1][2];
                    shell_symbol_A = "s";
                    } else if (i == 2){
                        exponent_A = c_exponents;
                        coefficient_A = c_2s_coefficients;
                        centerA_arma(0) = C_coords[0][0];
                        centerA_arma(1) = C_coords[0][1];
                        centerA_arma(2) = C_coords[0][2];
                        shell_symbol_A = "s";
                        
                        } else if (i == 3){
                            exponent_A = c_exponents;
                            coefficient_A = c_2p_coefficients;
                            centerA_arma(0) = C_coords[0][0];
                            centerA_arma(1) = C_coords[0][1];
                            centerA_arma(2) = C_coords[0][2];
                            shell_symbol_A = "p";
                            a = 0; // px
                            } else if (i == 4){
                            exponent_A = c_exponents;
                            coefficient_A = c_2p_coefficients;
                            centerA_arma(0) = C_coords[0][0];
                            centerA_arma(1) = C_coords[0][1];
                            centerA_arma(2) = C_coords[0][2];
                            shell_symbol_A = "p";
                            a = 1; // py
                            } else if (i == 5){
                            exponent_A = c_exponents;
                            coefficient_A = c_2p_coefficients;
                            centerA_arma(0) = C_coords[0][0];
                            centerA_arma(1) = C_coords[0][1];
                            centerA_arma(2) = C_coords[0][2];
                            shell_symbol_A = "p";
                            a = 2; // pz
                            } else if (i == 6) {
                                exponent_A = c_exponents;
                                coefficient_A = c_2s_coefficients;
                                centerA_arma(0) = C_coords[1][0];
                                centerA_arma(1) = C_coords[1][1];
                                centerA_arma(2) = C_coords[1][2];
                                shell_symbol_A = "s";
                            
                                } else if (i == 7) {
                                    exponent_A = c_exponents;
                                    coefficient_A = c_2p_coefficients;
                                    centerA_arma(0) = C_coords[1][0];
                                    centerA_arma(1) = C_coords[1][1];
                                    centerA_arma(2) = C_coords[1][2];
                                    shell_symbol_A = "p";
                                    a = 0; // px
                                    } else if (i == 8) {
                                    exponent_A = c_exponents;
                                    coefficient_A = c_2p_coefficients;
                                    centerA_arma(0) = C_coords[1][0];
                                    centerA_arma(1) = C_coords[1][1];
                                    centerA_arma(2) = C_coords[1][2];
                                    shell_symbol_A = "p";
                                    a = 1;
                                    } else if (i == 9) {
                                    exponent_A = c_exponents;
                                    coefficient_A = c_2p_coefficients;
                                    centerA_arma(0) = C_coords[1][0];
                                    centerA_arma(1) = C_coords[1][1];
                                    centerA_arma(2) = C_coords[1][2];
                                    shell_symbol_A = "p";
                                    a = 2;
                                    } else if (i == 10){
                                        exponent_A = h_exponents;
                                        coefficient_A = h_coefficients;
                                        centerA_arma(0) = H_coords[2][0];
                                        centerA_arma(1) = H_coords[2][1];
                                        centerA_arma(2) = H_coords[2][2];
                                        shell_symbol_A = "s";
                                        }
                            
            if (j == 1) {
                exponent_B = h_exponents;
                coefficient_B = h_coefficients;
                centerB_arma(0) = H_coords[1][0];
                centerB_arma(1) = H_coords[1][1];
                centerB_arma(2) = H_coords[1][2];
                shell_symbol_B = "s";

                } else if (j == 2){
                    exponent_B = c_exponents;
                    coefficient_B = c_2s_coefficients;
                    centerB_arma(0) = C_coords[0][0];
                    centerB_arma(1) = C_coords[0][1];
                    centerB_arma(2) = C_coords[0][2];
                    shell_symbol_B = "s";

                    } else if (j == 3) {
                        exponent_B = c_exponents;
                        coefficient_B = c_2p_coefficients;
                        centerB_arma(0) = C_coords[0][0];
                        centerB_arma(1) = C_coords[0][1];
                        centerB_arma(2) = C_coords[0][2];
                        shell_symbol_B = "p";
                        b = 0;
                        } else if (j == 4) {
                        exponent_B = c_exponents;
                        coefficient_B = c_2p_coefficients;
                        centerB_arma(0) = C_coords[0][0];
                        centerB_arma(1) = C_coords[0][1];
                        centerB_arma(2) = C_coords[0][2];
                        shell_symbol_B = "p";
                        b = 1;
                        } else if (j == 5) {
                        exponent_B = c_exponents;
                        coefficient_B = c_2p_coefficients;
                        centerB_arma(0) = C_coords[0][0];
                        centerB_arma(1) = C_coords[0][1];
                        centerB_arma(2) = C_coords[0][2];
                        shell_symbol_B = "p";
                        b = 2;
                        } else if (j == 6) {
                            exponent_B = c_exponents;
                            coefficient_B = c_2s_coefficients;
                            centerB_arma(0) = C_coords[1][0];
                            centerB_arma(1) = C_coords[1][1];
                            centerB_arma(2) = C_coords[1][2];
                            shell_symbol_B = "s";

                            } else if (j == 7) {
                                exponent_B = c_exponents;
                                coefficient_B = c_2p_coefficients;
                                centerB_arma(0) = C_coords[1][0];
                                centerB_arma(1) = C_coords[1][1];
                                centerB_arma(2) = C_coords[1][2];
                                shell_symbol_B = "p";
                                b = 0;
                                }  else if (j == 8) {
                                exponent_B = c_exponents;
                                coefficient_B = c_2p_coefficients;
                                centerB_arma(0) = C_coords[1][0];
                                centerB_arma(1) = C_coords[1][1];
                                centerB_arma(2) = C_coords[1][2];
                                shell_symbol_B = "p";
                                b = 1;
                                } else if (j == 9) {
                                exponent_B = c_exponents;
                                coefficient_B = c_2p_coefficients;
                                centerB_arma(0) = C_coords[1][0];
                                centerB_arma(1) = C_coords[1][1];
                                centerB_arma(2) = C_coords[1][2];
                                shell_symbol_B = "p";
                                b = 2;
                                } else if (j == 10) {
                                    exponent_B = h_exponents;
                                    coefficient_B = h_coefficients;
                                    centerB_arma(0) = H_coords[2][0];
                                    centerB_arma(1) = H_coords[2][1];
                                    centerB_arma(2) = H_coords[2][2];
                                    shell_symbol_B = "s";

                                    } else if (j == 11) {
                                        exponent_B = h_exponents;
                                        coefficient_B = h_coefficients;
                                        centerB_arma(0) = H_coords[3][0];
                                        centerB_arma(1) = H_coords[3][1];
                                        centerB_arma(2) = H_coords[3][2];
                                        shell_symbol_B = "s";

                                        }

            // within the two outer loops, there are two additional loops of k and l, each loop has 3 elements
            // general equationto calculate the basis functions: ∑∑h_coeff * c_coeff * norm_const * norm_const * overlap_integral
            // overlap_integral calculation depends on Shell definition
            double basis_function_ij = 0.0;
            for (int k=0; k < 3; k++){
                for (int l=0; l < 3; l++){

                    Shell shellA(centerA_arma, shell_symbol_A, exponent_A[k]);
                    Shell shellB(centerB_arma, shell_symbol_B, exponent_B[l]);

                    // Get the quantum numbers of the two shells
                    arma::imat quantumNumberA = shellA.getQuantumNumber();
                    arma::imat quantumNumberB = shellB.getQuantumNumber();

                    // double overlap = calculateOverlapIntegral_3D(centerA_arma, centerB_arma, exponent_A, exponent_B, quantumNumberA(k, l), quantumNumberB(k, l), k, l);
                    
                    if (shell_symbol_A == "s" && shell_symbol_B == "s") {
                        double overlap = calculateOverlapIntegral_ss(shellA, shellB);
                        // cout << "overlap: " << overlap << endl;
                        basis_function_ij += coefficient_A[k] * coefficient_B[l] * norm_const[i][k] * norm_const[i][l] * overlap;
                    } else if (shell_symbol_A == "p" && shell_symbol_B == "p") {
                        arma::mat overlap = calculateOverlapIntegral_pp(shellA, shellB);
                        // cout << "overlap: " << overlap(k, k) << endl;
                        basis_function_ij += coefficient_A[k] * coefficient_B[l] * norm_const[i][k] * norm_const[i][l] * overlap(a, b);
                    } else if (shell_symbol_A == "s" && shell_symbol_B == "p") {
                        std::vector<double> overlap = calculateOverlapIntegral_sp(shellA, shellB);
                        // cout << "overlap: " << overlap[k] << endl;
                        // cout << "overlap: " << overlap[l] << endl;
                        basis_function_ij += coefficient_A[k] * coefficient_B[l] * norm_const[i][k] * norm_const[i][l] * overlap[b];
                    } else if (shell_symbol_A == "p" && shell_symbol_B == "s") {
                        std::vector<double> overlap = calculateOverlapIntegral_sp(shellA, shellB);
                        // cout << "overlap: " << overlap[k] << endl;
                        // cout << "overlap: " << overlap[l] << endl;
                        basis_function_ij += coefficient_A[k] * coefficient_B[l] * norm_const[i][k] * norm_const[i][l] * overlap[a];
                    }
                }
            }

        // fill in matrix elements of basis_function_ij
        basis_functions(i, j) = basis_function_ij;
        basis_functions(j, i) = basis_function_ij;
        }
    }

    // define when i = j, basis function is 1
    for (int i=0; i<12; i++){
        basis_functions(i, i) = 1.0;
    }

    return basis_functions;
}

// calculate the hamiltonian of C2H4 as a 12*12 matrix
    // ({"H_s", -13.6});
    // ({"C_s", -21.4});
    // ({"C_px", -11.4});
    // ({"C_py", -11.4});
    // ({"C_pz", -11.4});
    // H(i, j) = K / 2 * (H_i + H_j) * S(i, j)
arma::mat molecule::calculate_hamiltonian_C2H4(arma::mat overlap_matrix_C2H4){
    
    // create a 12*12 matrix to store the hamiltonians
    arma::mat hamiltonian = arma::zeros<arma::mat>(12, 12);
    double K = 1.75;
    double h_1;
    double h_2;

    for (int i=0; i<11; i++){
        for (int j = i+1; j < 12; j++){
            if (i == 0){
                h_1 = -13.6;
                    
                } else if (i == 1) {
                    h_1 = -13.6;

                    } else if (i == 2){
                        h_1 = -21.4;
                    
                        } else if (i == 3 || i == 4 || i == 5){
                            h_1 = -11.4;
                        
                            } else if (i == 6) {
                                h_1 = -21.4;
                            
                                } else if (i == 7 || i == 8 || i == 9) {
                                    h_1 = -11.4;
                                    
                                    } else if (i == 10){
                                        h_1 = -13.6;
                                        }
                            
            if (j == 1) {
                h_2 = -13.6;

                } else if (j == 2){
                    h_2 = -21.4;

                    } else if (j == 3 || j == 4 || j == 5) {
                        h_2 = -11.4;

                        } else if (j == 6) {
                            h_2 = -21.4;

                            } else if (j == 7 || j == 8 || j == 9) {
                                h_2 = -11.4;

                                } else if (j == 10) {
                                    h_2 = -13.6;

                                    } else if (j == 11) {
                                        h_2 = -13.6;

                                        }

            // fill in matrix elements of basis_function_ij
            hamiltonian(i, j) = K / 2 * (h_1 + h_2) * overlap_matrix_C2H4(i, j);
            hamiltonian(j, i) = K / 2 * (h_1 + h_2) * overlap_matrix_C2H4(i, j);
            }
        }

        // define when i = j, basis function is 1
        for (int i=0; i<12; i++){
            hamiltonian(i, i) = h_1;
        }

        return hamiltonian;
}