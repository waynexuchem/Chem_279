// header files
# pragma once

#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <armadillo>
#include <string>

using namespace std;

class Shell {

public:
    
    arma::vec center;
    arma::imat quantumNumber;
    double exponent;
    
    // Constructor to initialize the shell with center R, shell type (s or p), and exponent (alpha)
    Shell(const arma::vec &center, const std::string& shellType, double exponent);

    // Function to get the quantum numbers
    arma::imat getQuantumNumber() const;

    // Function to get the center of the shell
    arma::vec getCenter() const;

    // Function to get the exponent alpha
    double getExponent() const;

};