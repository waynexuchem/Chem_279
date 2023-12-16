/*
Numerically evaluate the indefinite 1D overlap integral between two Gaussian type functions, GA(x) and GB(x).
GA(x) = (x-XA)^lA * exp[-alpha * (x-XA)^2]
GB(x) = (x-XB)^lB * exp[-beta * (x-XB)^2]
SxAB = ∫−∞->∞GA(x)⋅GB(x)dx

Integration method is the extended trapezoidal rule with equally spaced points.
*/

// import packages
#include "Utils.h"

#include <iostream>
#include <cmath>
#include <fstream> // for reading/writing files

using namespace std;

// Two s-type functions (l A = lB = 0) both centered at the origin, with an integration range of -5 to 5
void test_case_1_1() {
    double XA = 0.0;
    double alpha = 1;
    int lA = 0;
    double XB = 0.0;
    double beta = 1;
    int lB = 0;
    double lower_bound = -5.0;
    double higher_bound = 5.0;
    int num_intervals = 10000; // Maximum number of intervals
    double result = SxAB(gaussian_func, lower_bound, higher_bound, XA, alpha, lA, XB, beta, lB, num_intervals);

    cout << "The indefinite overlap integral between two s-type functions (l A = lB = 0) both centered at the origin, with an integration range of -5 to 5 is: " << result << endl;
}

// Two s-type functions (l A = lB = 0) both centered at the origin, with an integration range of -1 to 1
void test_case_1_2() {
    double XA = 0.0;
    double alpha = 1;
    int lA = 0;
    double XB = 0.0;
    double beta = 1;
    int lB = 0;
    double lower_bound = -1.0;
    double higher_bound = 1.0;
    int num_intervals = 10000; // Maximum number of intervals
    double result = SxAB(gaussian_func, lower_bound, higher_bound, XA, alpha, lA, XB, beta, lB, num_intervals);

    cout << "The indefinite overlap integral between two s-type functions (l A = lB = 0) both centered at the origin, with an integration range of -1 to 1 is: " << result << endl;
}

// Two s-type functions (l A = lB = 0) with an offset of 1.0, with an integration range of -5 to 5
void test_case_2() {
    double XA = 0.0;
    double alpha = 1;
    int lA = 0;
    double XB = 1.0;
    double beta = 1;
    int lB = 0;
    double lower_bound = -5.0;
    double higher_bound = 5.0;
    int num_intervals = 10000; // Maximum number of intervals

    vector<int> h_values;
    vector<double> integration_results;
    vector<double> estimated_errors;
    double result = SxAB(gaussian_func, lower_bound, higher_bound, XA, alpha, lA, XB, beta, lB, num_intervals);

    cout << "The indefinite overlap integral between two s-type functions (l A = lB = 0) with an offset of 1.0, with an integration range of -5 to 5 is: " << result << endl;
}

// An s-type function and a p-type function (lA = 0, lB = 1) both centered at the origin, with an integration range of -5 to 5
void test_case_3() {
    double XA = 0.0;
    double alpha = 1;
    int lA = 0;
    double XB = 0.0;
    double beta = 1;
    int lB = 1;
    double lower_bound = -5.0;
    double higher_bound = 5.0;
    int num_intervals = 10000; // Maximum number of intervals
    double result = SxAB(gaussian_func, lower_bound, higher_bound, XA, alpha, lA, XB, beta, lB, num_intervals);

    cout << "The indefinite overlap integral between an s-type function and a p-type function (lA = 0, lB = 1) both centered at the origin, with an integration range of -5 to 5 is: " << result << endl;
}

// An s-type function and a p-type function (lA = 0, lB = 1) with an offset of 1.0, with an integration range of -5 to 5
// Rebecca got -0.6
void test_case_4() {
    double XA = 0.0;
    double alpha = 1.0;
    int lA = 0;
    double XB = 1.0;
    double beta = 1.0;
    int lB = 1;
    double lower_bound = -5.0;
    double higher_bound = 5.0;
    int num_intervals = 10000; // Maximum number of intervals
    double result = SxAB(gaussian_func, lower_bound, higher_bound, XA, alpha, lA, XB, beta, lB, num_intervals);

    cout << "The indefinite overlap integral between an s-type function and a p-type function (lA = 0, lB = 1) with an offset of 1.0, with an integration range of -5 to 5 is: " << result << endl;
}



int main() {
    test_case_1_1();
    test_case_1_2();
    test_case_2();
    test_case_3();
    test_case_4();
    return 0;
}