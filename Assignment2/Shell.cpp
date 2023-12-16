#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <armadillo>
#include <string>

#include "Shell.h" // Include the Shell class header

Shell::Shell(const arma::vec & R, const std::string& shellType, double alpha) : center(R), exponent(alpha) {
    // Determine quantum number based on shellType
    if (shellType == "s") {
        quantumNumber = arma::zeros<arma::imat>(3, 3); // matrix of [0, 0, 0; 0, 0, 0; 0, 0, 0]
    } else if (shellType == "p") {
        quantumNumber = arma::eye<arma::imat>(3, 3);  // matrix of [1, 0, 0; 0, 1, 0; 0, 0, 1]
    }
}

arma::imat Shell::getQuantumNumber() const {
    return quantumNumber;
}

arma::vec Shell::getCenter() const
{
    return center;
}

double Shell::getExponent() const {
    return exponent;
}