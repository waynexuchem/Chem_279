// import packages
#include "Shell.h"
#include "Utils.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

// Define Gaussian_type_function(x)
double gaussian_func(double x, double center, double exponent, int power) {
    return pow(x - center, power) * exp(-exponent * pow(x - center, 2));
}

// Define the indefinite overlap integral
// SxAB = ∫−∞->∞GA(x)⋅GB(x)dx
// Name the funtion trapezoidal_integral
double SxAB(double (*func) (double, double, double, int), double lower_bound, double higher_bound, double XA, double alpha, int lA, double XB, double beta, int lB, int num_intervals) {
    double interval = (higher_bound - lower_bound) / num_intervals;
    double result = 0.0;

    for (int i = 0; i <= num_intervals; ++i) {
        double x = lower_bound + i * interval;
        double value = func(x, XA, alpha, lA) * func(x, XB, beta, lB);

        if (i == 0 || i == num_intervals) {
            result += value;
        } else {
            result += 2 * value;
        }
    }
    return 0.5 * interval * result;
}

// // Define central difference approximation of the second derivative of SxAB
// // f'' = [f(x+h) - 2f(x) + f(x-h)] / h^2
// double d2SxAB(double (*func) (double, double, double, int), double x, double h, double XA, double alpha, int lA, double XB, double beta, int lB) {
//     double f_plus = func(x + h, XA, alpha, lA) * func(x + h, XB, beta, lB);
//     double f_minus = func(x - h, XA, alpha, lA) * func(x - h, XB, beta, lB);
//     double f = func(x, XA, alpha, lA) * func(x, XB, beta, lB);

//     return (f_plus - 2 * f + f_minus) / pow(h, 2);
// }

// // Define write_convergence_interval function to write integration results and error to a file
// void write_convergence_interval(const std::string& file_path, std::vector<int> h_values, 
//                             std::vector<double> integration_results, std::vector<double> estimated_errors){

//         // Opens up a file stream for output
//         std::ofstream outfile(file_path);

//         outfile.clear();             // Clear the EOF flag

//         // Check that it was successfully opened
//         if(!outfile.is_open())
//         {   
//             throw std::runtime_error("File path in write_truncation_error does not exist!");
//         }

//         for (int i = 0; i < h_values.size(); i++){
//             outfile << h_values[i] << "\t" << integration_results[i] << "\t" << estimated_errors[i] << std::endl;
//         }

//     }

// Write a function to calculate factorial
int factorial(int n) {
    if (n == 0) {
        return 1;
    } else {
        return n * factorial(n-1);
    }
}

// Write a function to calculate double factorial
// n!! = n(n-2)(n-4)...(4)(2) if n is even
// n!! = n(n-2)(n-4)...(3)(1) if n is odd
// n!! = (n+2)!!/(n+2) if n < 0 and n is odd
int doubleFactorial(int n) {
    if (n == 0 || n == 1) {
        return 1;
    }

    if (n < 0 && n % 2 != 0) {
        return doubleFactorial(n + 2) / (n + 2);
    }

    int result = 1;
    for (int i = n; i >= 1; i -= 2) {
        result *= i;
    }
    return result;
}


// Write a function to calculate integer bionomial_coefficient (Equation 2.8)
// m choose n = m! / [n! * (m - n)!]
// bionomial_coefficient = sup! / [sub! * (sup - sub)!]
int bionomial_coefficient(int m, int n) {
    return factorial(m) / (factorial(n) * factorial(m - n));
}

// Write a function to calculate center of the product shell
// R_P = (α*R_A + β*R_B) / (α+β)
arma::vec calculateCenter(arma::vec centerA, arma::vec centerB, double alpha, double beta) {
    
    arma::vec centerP(3);
    for (int i = 0; i < 3; ++i) {
        centerP[i] = (alpha * centerA[i] + beta * centerB[i]) / (alpha + beta);
    }

    return centerP;
}

// Write a function to calculate the analytical overlap integral from all x, y, z components between two s-shells
// with given quantum numbers of the two shells
// S^AB_x = exp [−α*β*(X_A −X_B )^2 / (α+β)] * √(π / (α+β)) * {∑^lA_i=0[∑^lB_j=0 (bionomial_coefficient(lA, i) * bionomial_coefficient(lB, j) * (i + j - 1)!! * (X_P - X_A)^(lA-i) * (X_P - X_B)^(lb-i) / [2 * (α+β)]^[(i+j)/2]]}
// The result is a scalar
double calculateOverlapIntegral_ss(Shell ShellA, Shell ShellB) {
    // Get the quantum numbers of the two shells
    arma::imat quantumNumberA = ShellA.getQuantumNumber();
    arma::imat quantumNumberB = ShellB.getQuantumNumber();

    // Get the center and exponent of the two shells
    arma::vec centerA = ShellA.getCenter();
    arma::vec centerB = ShellB.getCenter();
    double alpha = ShellA.getExponent();
    double beta = ShellB.getExponent();

    // Calculate the center of the product shell
    arma::vec centerP = calculateCenter(centerA, centerB, alpha, beta);

    // Define the quantum numbers of the two shells
    int lA = quantumNumberA(0, 0);
    int lB = quantumNumberB(0, 0);

    // Calculate the overlap integral
    double result = 0.0;

    // Calculate component 1 for x, y, z: exp [−α*β*(X_A −X_B )^2 / (α+β)]
    double component1_x = exp(-alpha * beta * pow(centerA[0] - centerB[0], 2) / (alpha + beta));
    double component1_y = exp(-alpha * beta * pow(centerA[1] - centerB[1], 2) / (alpha + beta));
    double component1_z = exp(-alpha * beta * pow(centerA[2] - centerB[2], 2) / (alpha + beta));

    // Calculate component 2 for x, y, z: √(π / (α+β))
    double component2 = sqrt(M_PI / (alpha + beta));

    // Calculate the double integral component for x, y, z
    double component3_x = 0.0;
    double component3_y = 0.0;
    double component3_z = 0.0;

    for (int i=0; i<=lA; i++){
        for (int j=0; j<=lB; j++){
            component3_x += bionomial_coefficient(lA, i) * bionomial_coefficient(lB, j) * doubleFactorial(i + j - 1) * pow(centerP[0] - centerA[0], lA - i) * pow(centerP[0] - centerB[0], lB - j) / pow(2 * (alpha + beta), (i + j) / 2);
            component3_y += bionomial_coefficient(lA, i) * bionomial_coefficient(lB, j) * doubleFactorial(i + j - 1) * pow(centerP[1] - centerA[1], lA - i) * pow(centerP[1] - centerB[1], lB - j) / pow(2 * (alpha + beta), (i + j) / 2);
            component3_z += bionomial_coefficient(lA, i) * bionomial_coefficient(lB, j) * doubleFactorial(i + j - 1) * pow(centerP[2] - centerA[2], lA - i) * pow(centerP[2] - centerB[2], lB - j) / pow(2 * (alpha + beta), (i + j) / 2);
        }
    }

    // Calculate the overlap integral for x, y, z
    double overlap_x = component1_x * component2 * component3_x;
    double overlap_y = component1_y * component2 * component3_y;
    double overlap_z = component1_z * component2 * component3_z;

    // Calculate the overlap integral
    result = overlap_x * overlap_y * overlap_z;

    return result;
}

// Write a function to calculate the analytical overlap integral from all x, y, z components between an s-shell and a p-shell
// with given quantum numbers of the two shells
// S^AB_x = exp [−α*β*(X_A −X_B )^2 / (α+β)] * √(π / (α+β)) * {∑^lA_i=0[∑^lB_j=0 (bionomial_coefficient(lA, i) * bionomial_coefficient(lB, j) * (i + j - 1)!! * (X_P - X_A)^(lA-i) * (X_P - X_B)^(lb-i) / [2 * (α+β)]^[(i+j)/2]]}
// The result is a vector of size 3
std::vector<double> calculateOverlapIntegral_sp(Shell sShell, Shell pShell) {
    
    // Get the quantum numbers of the two shells
    arma::imat quantumNumberS = sShell.getQuantumNumber();
    arma::imat quantumNumberP = pShell.getQuantumNumber();

    // Get the center and exponent of the two shells
    arma::vec centerA = sShell.getCenter();
    arma::vec centerB = pShell.getCenter();
    double alpha = sShell.getExponent();
    double beta = pShell.getExponent();

    // Calculate the center of the product shell
    arma::vec centerP = calculateCenter(centerA, centerB, alpha, beta);

    // Calculate component 1 for x, y, z: exp [−α*β*(X_A −X_B )^2 / (α+β)]
    double component1_x = exp(-alpha * beta * pow(centerA[0] - centerB[0], 2) / (alpha + beta));
    double component1_y = exp(-alpha * beta * pow(centerA[1] - centerB[1], 2) / (alpha + beta));
    double component1_z = exp(-alpha * beta * pow(centerA[2] - centerB[2], 2) / (alpha + beta));

    // Calculate component 2 for x, y, z: √(π / (α+β))
    double component2 = sqrt(M_PI / (alpha + beta));

    std::vector<double> result(3, 0.0);

    // Define the quantum numbers of the two shells
    int lA = quantumNumberS(0, 0);
    arma::vec lBx = arma::conv_to<arma::vec>::from(quantumNumberP.row(0)); // Extract the first row [1, 0, 0] corresponding to l, m, n
    arma::vec lBy = arma::conv_to<arma::vec>::from(quantumNumberP.row(1)); // Extract the second row [0, 1, 0] corresponding to l, m, n
    arma::vec lBz = arma::conv_to<arma::vec>::from(quantumNumberP.row(2)); // Extract the third row [0, 0, 1] corresponding to l, m, n

    // Calculate the double integral component for x, y, z axes (denoted by l = 0, 1, 2) 
    // with quantum numbers lA & lB
    for (int l = 0; l < 3; ++l) {
        std::vector<double> component3(3, 1.0); // Initialize component3 vector for x, y, or z to 1.0

        arma::vec lB; // Initialize lB to the appropriate quantum numbers based on the current orbital
        if (l == 0) {
            lB = lBx;
        } else if (l == 1) {
            lB = lBy;
        } else if (l == 2) {
            lB = lBz;
        }

        // Another loop to go through all the quantum numbers of the p-shell within lB
        // Calculate sp_l, sp_m, sp_n
        for (int q = 0; q < 3; q++){
            double sum_factor = 0.0;
            for (int i = 0; i <= lA; ++i) {
                    for (int j = 0; j <= lB(q); ++j) {
                        if ((i + j) % 2 != 0) {
                            continue;
                        }

                        sum_factor += bionomial_coefficient(lA, i) * bionomial_coefficient(lB(q), j) * doubleFactorial(i + j - 1) *
                            pow(centerP[l] - centerA[l], lA - i) * pow(centerP[l] - centerB[l], lB(q) - j) /
                            pow(2 * (alpha + beta), (i + j) / 2);
                        // cout << "sum_factor: " << sum_factor << "\n";

                    }
                }
            component3[l] *= sum_factor;
            // cout << "component3[("<<l<<")]: " << component3[l] << "\n";
            }
        

        // Calculate the overlap integral for x, y, or z
        if (l == 0) {
            result[l] = component1_x * component2 * component3[l];
        } else if (l == 1) {
            result[l] = component1_y * component2 * component3[l];
        } else if (l == 2) {
            result[l] = component1_z * component2 * component3[l];
        }
        
    }

    return result;
}

// Write a function to calculate the analytical overlap integral from all x, y, z components between two p-shells
// with given quantum numbers of the two shells
// S^AB_x = exp [−α*β*(X_A −X_B )^2 / (α+β)] * √(π / (α+β)) * {∑^lA_i=0[∑^lB_j=0 (bionomial_coefficient(lA, i) * bionomial_coefficient(lB, j) * (i + j - 1)!! * (X_P - X_A)^(lA-i) * (X_P - X_B)^(lb-i) / [2 * (α+β)]^[(i+j)/2]]}
// The result is a matrix of [3, 3]
// Shape is:
// [[xx, xy, xz],
//  [yx, yy, yz],
//  [zx, zy, zz]]
arma::mat calculateOverlapIntegral_pp(Shell ShellA, Shell ShellB) {
    // Input shellA is the p-shell
    // Input shellB is the p-shell
    
    // Get the quantum numbers of the two shells
    arma::imat quantumNumberA = ShellA.getQuantumNumber();
    arma::imat quantumNumberB = ShellB.getQuantumNumber();

    // Get the center and exponent of the two shells
    arma::vec centerA = ShellA.getCenter();
    arma::vec centerB = ShellB.getCenter();
    double alpha = ShellA.getExponent();
    double beta = ShellB.getExponent();

    // Calculate the center of the product shell
    arma::vec centerP = calculateCenter(centerA, centerB, alpha, beta);

    // Define the quantum numbers of the two shells
    arma::vec lAx = arma::conv_to<arma::vec>::from(quantumNumberA.row(0)); // Extract the first row [1, 0, 0] corresponding to l, m, n
    arma::vec lAy = arma::conv_to<arma::vec>::from(quantumNumberA.row(1)); // Extract the second row [0, 1, 0] corresponding to l, m, n
    arma::vec lAz = arma::conv_to<arma::vec>::from(quantumNumberA.row(2)); // Extract the third row [0, 0, 1] corresponding to l, m, n
    arma::vec lBx = arma::conv_to<arma::vec>::from(quantumNumberB.row(0)); // Extract the first row [1, 0, 0] corresponding to l, m, n
    arma::vec lBy = arma::conv_to<arma::vec>::from(quantumNumberB.row(1)); // Extract the second row [0, 1, 0] corresponding to l, m, n
    arma::vec lBz = arma::conv_to<arma::vec>::from(quantumNumberB.row(2)); // Extract the third row [0, 0, 1] corresponding to l, m, n

    // Calculate the overlap integral
    arma::mat result(3, 3);

    // Calculate component 1 for x, y, z: exp [−α*β*(X_A −X_B )^2 / (α+β)]
    double component1_x = exp(-alpha * beta * pow(centerA[0] - centerB[0], 2) / (alpha + beta));
    double component1_y = exp(-alpha * beta * pow(centerA[1] - centerB[1], 2) / (alpha + beta));
    double component1_z = exp(-alpha * beta * pow(centerA[2] - centerB[2], 2) / (alpha + beta));

    // Calculate component 2 for x, y, z: √(π / (α+β))
    double component2 = sqrt(M_PI / (alpha + beta));

    // Calculate the double integral component for x, y, z axes with quantum numbers lA & lB
    for (int l=0; l<3; l++){
        
        arma::vec lA; // Initialize lA to the appropriate quantum numbers based on the current orbital

        if (l == 0) {
            lA = lAx;
        } else if (l == 1) {
            lA = lAy;
        } else if (l == 2) {
            lA = lAz;
        }
        
        

        for (int m=0; m<3; m++){

            arma::vec lB; // Initialize lA to the appropriate quantum numbers based on the current orbital

            if (m == 0) {
                lB = lBx;
            } else if (m == 1) {
                lB = lBy;
            } else if (m == 2) {
                lB = lBz;
            }

        double component3 = 1.0; // Initialize component3 to 1.0
        // Another loop to go through all the quantum numbers of the p-shell within lA and lB
        // Calculate pp_ll, pp_mm, pp_nn
        for (int q = 0; q < 3; q++){
            double sum_factor = 0.0;
            for (int i = 0; i <= lA(q); ++i) {
                for (int j = 0; j <= lB(q); ++j) {
                    if ((i + j) % 2 != 0) {
                        continue;
                    }

                    sum_factor += bionomial_coefficient(lA(q), i) * bionomial_coefficient(lB(q), j) * doubleFactorial(i + j - 1) *
                        pow(centerP[l] - centerA[l], lA(q) - i) * pow(centerP[l] - centerB[l], lB(q) - j) /
                        pow(2 * (alpha + beta), (i + j) / 2);

                }
            }
            component3 *= sum_factor;
            // cout << "component3["<<l<<", "<<m<<"]: " << component3 << "\n";
        }

        // Calculate the overlap integral and append to the result vector
        // Calculate the overlap integral for x, y, or z
        if (m == 0) {
            result(l,m) = component1_x * component2 * component3;
        } else if (m == 1) {
            result(l,m) = component1_y * component2 * component3;
        } else if (m == 2) {
            result(l,m) = component1_z * component2 * component3;
        }
        }
    }

    return result;
}