// header files
# pragma once
#include "Shell.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <armadillo>
#include <string>

using namespace std;

// Define Gaussian_type_function(x)
double gaussian_func(double x, double center, double exponent, int power);

// Define the indefinite overlap integral
// SxAB = ∫−∞->∞GA(x)⋅GB(x)dx
double SxAB(double (*func) (double, double, double, int), double lower_bound, double higher_bound, double XA, double alpha, int lA, double XB, double beta, int lB, int num_intervals);

// Define central difference approximation of the second derivative of SxAB
// f'' = [f(x+h) - 2f(x) + f(x-h)] / h^2
double d2SxAB(double (*func) (double, double, double, int), double x, double h, double XA, double alpha, int lA, double XB, double beta, int lB);

// Define write_convergence_interval function to write integration results and error to a file
void write_convergence_interval(const std::string& file_path, std::vector<int> h_values, 
                            std::vector<double> integration_results, std::vector<double> estimated_errors);

// Write a function to calculate factorial
int factorial(int n);

// Write a function to calculate double factorial
// n!! = n(n-2)(n-4)...(4)(2) if n is even
// n!! = n(n-2)(n-4)...(3)(1) if n is odd
// n!! = (n+2)!!/(n+2) if n < 0 and n is odd
int doubleFactorial(int n);

// Write a function to calculate integer prefactor
// prefactor = sup! / [sub! * (sup - sub)!]
int prefactor(int superscript, int subscript);

// Write a function to calculate center of the product shell
// R_P = (α*R_A + β*R_B) / (α+β)
std::vector<double> calculateCenter(std::vector<double> centerA, std::vector<double> centerB, double alpha, double beta);

// Write a function to calculate the analytical overlap integral from all x, y, z components between two s-shells
// with given quantum numbers of the two shells
// S^AB_x = exp [−α*β*(X_A −X_B )^2 / (α+β)] * √(π / (α+β)) * {∑^lA_i=0[∑^lB_j=0 (prefactor(lA, i) * prefactor(lB, j) * (i + j - 1)!! * (X_P - X_A)^(lA-i) * (X_P - X_B)^(lb-i) / [2 * (α+β)]^[(i+j)/2]]}
// The result is a scalar
double calculateOverlapIntegral_ss(Shell ShellA, Shell ShellB);

// Write a function to calculate the analytical overlap integral from all x, y, z components between an s-shell and a p-shell
// with given quantum numbers of the two shells
// S^AB_x = exp [−α*β*(X_A −X_B )^2 / (α+β)] * √(π / (α+β)) * {∑^lA_i=0[∑^lB_j=0 (prefactor(lA, i) * prefactor(lB, j) * (i + j - 1)!! * (X_P - X_A)^(lA-i) * (X_P - X_B)^(lb-i) / [2 * (α+β)]^[(i+j)/2]]}
// The result is a vector of size 3
std::vector<double> calculateOverlapIntegral_sp(Shell ShellA, Shell ShellB);

// Write a function to calculate the analytical overlap integral from all x, y, z components between two p-shells
// with given quantum numbers of the two shells
// S^AB_x = exp [−α*β*(X_A −X_B )^2 / (α+β)] * √(π / (α+β)) * {∑^lA_i=0[∑^lB_j=0 (prefactor(lA, i) * prefactor(lB, j) * (i + j - 1)!! * (X_P - X_A)^(lA-i) * (X_P - X_B)^(lb-i) / [2 * (α+β)]^[(i+j)/2]]}
// The result is a matrix of [3, 3]
// Shape is:
// [[xx, xy, xz],
//  [yx, yy, yz],
//  [zx, zy, zz]]
arma::mat calculateOverlapIntegral_pp(Shell ShellA, Shell ShellB);

// Split the product of integral into 3 dimensions
// and calculate each dimension separately
double calculateOverlapIntegral_1D(int dimension, arma::vec centerA, arma::vec centerB, double alpha, double beta, int lA, int lB);