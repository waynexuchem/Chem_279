#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <armadillo>
#include <string>

#include "Shell.h" // Include the Shell class header
#include "Utils.h"  // import Utils package
#include "molecule.h"

void test_case_1() {
    molecule H2("/Users/wangxu/Documents/Chem_279/Assignment/Assignment3/Rebecca_Ref/H2.txt");
    Coordinates H_coords = H2.H_coords;
    Coordinates C_coords = H2.C_coords;

    // print number of basis functions
    std::cout << "Number of basis functions: " << H2.calculate_num_basis_functions(C_coords, H_coords) << std::endl;

    // print center of the two H atoms
    const AtomCoord centerH1 = H2.calculateCenter(H_coords[0], H_coords[1], h_exponents[0], h_exponents[0]);
    std::cout << "Center of the two H atoms: (" << centerH1[0] << ", " << centerH1[1] << ", " << centerH1[2] << ")" << std::endl;

    // print normalization constants for the two H atoms (S)
    const auto [normH1, normH2, normH3] = H2.normalization_constant_S(h_exponents, centerH1);
    std::cout << "Normalization constants for the two H atoms (S): " << std::endl;
    std::cout << "normH1 = " << normH1 << std::endl;
    std::cout << "normH2 = " << normH2 << std::endl;
    std::cout << "normH3 = " << normH3 << std::endl;

    // print off-diagonal elements of the overlap matrix (S)
    std::cout << "Off-diagonal elements of the overlap matrix (S): " << std::endl;
    std::cout << "S^HH_1 = " << H2.calculate_basis_functions_S(h_exponents, h_coefficients, H_coords[0], H_coords[1]) << std::endl;

    // Build an armadillo matrix for the overlap matrix (S)
    int mat_dim = H2.calculate_num_basis_functions(C_coords, H_coords);
    arma::mat S = arma::zeros<arma::mat>(mat_dim, mat_dim);
    S(0, 0) = 1.0;
    S(0, 1) = H2.calculate_basis_functions_S(h_exponents, h_coefficients, H_coords[0], H_coords[1]);
    S(1, 0) = S(0, 1);
    S(1, 1) = 1.0;
    std::cout << "Overlap matrix (S): " << std::endl;
    S.print();

    // build an armadillo matrix for Hamiltonian matrix (H)
    arma::mat H = arma::zeros<arma::mat>(mat_dim, mat_dim);
    H(0, 0) = -13.6;
    H(0, 1) = H2.calculate_hamiltonian_S(h_exponents, h_coefficients, H_coords[0], H_coords[1]);
    H(1, 0) = H(0, 1);
    H(1, 1) = -13.6;
    std::cout << "Hamiltonian matrix (H): " << std::endl;
    H.print();

    // diagonalize S(overlap) matrix -> gives us eigenvectors (matrix U) and eigenvalues (vector s)
    arma::vec s;
    arma::mat U;

    arma::eig_sym(s, U, S);  // eigen decomposition of S as a symmetric matrix
    s.print("Eigenvalues of S: ");
    U.print("Eigenvectors of S: ");

    // Get inverse square root of  of eigen values and put them in a diagonal matrix
    for (int i = 0; i < s.size(); i++){
        s(i) = 1.0 / sqrt(s(i));
    }
    arma::mat s_diag_mat = arma::diagmat(s);
    arma::mat X = U * s_diag_mat * U.t(); // X = U * s^(-1/2) * U^T
    X.print("X matrix: ");

    // Transform H matrix to the new basis
    arma::mat H_prime = X.t() * H * X;
    H_prime.print("H_prime matrix: ");

    // Generate the molecular orbital coefficients C_mat: C =XV
    arma::mat C_mat = X * U;
    C_mat.print("C_mat matrix: ");

    // Calculate the energy eigenvector of the system: HV =Vε
    arma::vec epsilon;
    arma::eig_sym(epsilon, U, H_prime);  // eigen decomposition of H_prime as a symmetric matrix
    epsilon.print("Eigenvalues of H_prime: ");

    // Calculate the energy of the system: E = ∑2εi
    double E = 2 * epsilon(0);

    std::cout << "Energy of the system: " << E << std::endl;

}

void test_case_2() {
    molecule C2H2("/Users/wangxu/Documents/Chem_279/Assignment/Assignment3/Rebecca_Ref/C2H2.txt");
    Coordinates H_coords = C2H2.H_coords;
    Coordinates C_coords = C2H2.C_coords;

    // print values of H_coords
    std::cout << "Values of H_coords: " << std::endl;
    for (int i = 0; i < H_coords.size(); i++) {
        for (int j = 0; j < H_coords[i].size(); j++) {
            std::cout << H_coords[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // print values of C_coords
    std::cout << "Values of C_coords: " << std::endl;
    for (int i = 0; i < C_coords.size(); i++) {
        for (int j = 0; j < C_coords[i].size(); j++) {
            std::cout << C_coords[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // print number of basis functions
    std::cout << "Number of basis functions: " << C2H2.calculate_num_basis_functions(C_coords, H_coords) << std::endl;

    // print normalization constant matrix of C2H2 as a 10*3 matrix
    std::cout << "Normalization constant matrix of C2H2: " << std::endl;
    std::vector<std::vector<double>> norm_const = C2H2.C2H2_norm_const(h_exponents, c_exponents, H_coords, C_coords);
    for (int i = 0; i < norm_const.size(); i++) {
        for (int j = 0; j < norm_const[i].size(); j++) {
            std::cout << norm_const[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // print basis function matrix of C2H2 as a 10*10 matrix
    arma::mat S = C2H2.calculate_basis_function_C2H2(h_exponents, c_exponents, h_coefficients, c_2s_coefficients, c_2p_coefficients, H_coords, C_coords);
    S.print("Overlap matrix of C2H2: ");

    // print Hamiltonian matrix of C2H2 as a 10*10 matrix
    arma::mat H = C2H2.calculate_hamiltonian_C2H2(S);
    H.print("Hamiltonian matrix of C2H2: ");

    // diagonalize S(overlap) matrix -> gives us eigenvectors (matrix U) and eigenvalues (vector s)
    arma::vec s;
    arma::mat U;

    arma::eig_sym(s, U, S);  // eigen decomposition of S as a symmetric matrix
    s.print("Eigenvalues of S: ");
    U.print("Eigenvectors of S: ");

    // Get inverse square root of  of eigen values and put them in a diagonal matrix
    for (int i = 0; i < s.size(); i++){
        s(i) = 1.0 / sqrt(s(i));
    }
    arma::mat s_diag_mat = arma::diagmat(s);
    arma::mat X = U * s_diag_mat * U.t(); // X = U * s^(-1/2) * U^T
    X.print("X matrix: ");

    // Transform H matrix to the new basis
    arma::mat H_prime = X.t() * H * X;
    H_prime.print("H_prime matrix: ");

    // Generate the molecular orbital coefficients C_mat: C =XV
    arma::mat C_mat = X * U;
    C_mat.print("C_mat matrix: ");

    // Calculate the energy eigenvector of the system: HV =Vε
    arma::vec epsilon;
    arma::eig_sym(epsilon, U, H_prime);  // eigen decomposition of H_prime as a symmetric matrix
    epsilon.print("Eigenvalues of H_prime: ");

    // Calculate the energy of the system: E = ∑2εi
    double E = 2 * epsilon(0);

    std::cout << "Energy of the system: " << E << std::endl;
}

void test_case_3() {
    molecule C2H4("/Users/wangxu/Documents/Chem_279/Assignment/Assignment3/Rebecca_Ref/C2H4.txt");
    Coordinates H_coords = C2H4.H_coords;
    Coordinates C_coords = C2H4.C_coords;

    // print values of H_coords
    std::cout << "Values of H_coords: " << std::endl;
    for (int i = 0; i < H_coords.size(); i++) {
        for (int j = 0; j < H_coords[i].size(); j++) {
            std::cout << H_coords[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // print values of C_coords
    std::cout << "Values of C_coords: " << std::endl;
    for (int i = 0; i < C_coords.size(); i++) {
        for (int j = 0; j < C_coords[i].size(); j++) {
            std::cout << C_coords[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // print number of basis functions
    std::cout << "Number of basis functions: " << C2H4.calculate_num_basis_functions(C_coords, H_coords) << std::endl;

    // print normalization constant matrix of C2H4 as a 12*3 matrix
    std::cout << "Normalization constant matrix of C2H4: " << std::endl;
    std::vector<std::vector<double>> norm_const = C2H4.C2H4_norm_const(h_exponents, c_exponents, H_coords, C_coords);
    for (int i = 0; i < norm_const.size(); i++) {
        for (int j = 0; j < norm_const[i].size(); j++) {
            std::cout << norm_const[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // print basis function matrix of C2H4 as a 12*12 matrix
    arma::mat S = C2H4.calculate_basis_function_C2H4(h_exponents, c_exponents, h_coefficients, c_2s_coefficients, c_2p_coefficients, H_coords, C_coords);
    S.print("Overlap matrix of C2H4: ");

    // print Hamiltonian matrix of C2H2 as a 10*10 matrix
    arma::mat H = C2H4.calculate_hamiltonian_C2H4(S);
    H.print("Hamiltonian matrix of C2H4: ");

    // diagonalize S(overlap) matrix -> gives us eigenvectors (matrix U) and eigenvalues (vector s)
    arma::vec s;
    arma::mat U;

    arma::eig_sym(s, U, S);  // eigen decomposition of S as a symmetric matrix
    s.print("Eigenvalues of S: ");
    U.print("Eigenvectors of S: ");

    // Get inverse square root of  of eigen values and put them in a diagonal matrix
    for (int i = 0; i < s.size(); i++){
        s(i) = 1.0 / sqrt(s(i));
    }
    arma::mat s_diag_mat = arma::diagmat(s);
    arma::mat X = U * s_diag_mat * U.t(); // X = U * s^(-1/2) * U^T
    X.print("X matrix: ");

    // Transform H matrix to the new basis
    arma::mat H_prime = X.t() * H * X;
    H_prime.print("H_prime matrix: ");

    // Generate the molecular orbital coefficients C_mat: C =XV
    arma::mat C_mat = X * U;
    C_mat.print("C_mat matrix: ");

    // Calculate the energy eigenvector of the system: HV =Vε
    arma::vec epsilon;
    arma::eig_sym(epsilon, U, H_prime);  // eigen decomposition of H_prime as a symmetric matrix
    epsilon.print("Eigenvalues of H_prime: ");

    // Calculate the energy of the system: E = ∑2εi
    double E = 2 * epsilon(0);

    std::cout << "Energy of the system: " << E << std::endl;
}

int main() {
    test_case_1();
    test_case_2();
    test_case_3();
    return 0;
}