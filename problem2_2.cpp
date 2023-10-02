#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <armadillo>
#include <string>

#include "Shell.h" // Include the Shell class header
#include "Utils.h"  // import Utils package


// Test case 1 to evaluate the overlap integral between two s-type functions (l A = lB = 0) both centered at the origin
void test_case_1() {
    std::vector<double> centerA = {0.0, 0.0, 0.0};
    std::vector<double> centerB = {0.0, 0.0, 0.0};
    double alpha = 1.0;
    double beta = 1.0;

    // Create Shell objects for s-shells A and B
    Shell shellA(centerA, "s", alpha);
    Shell shellB(centerB, "s", beta);

    // Calculate the overlap integral
    double overlap = calculateOverlapIntegral_ss(shellA, shellB);

    std::cout << "Overlap integral between two s-shells is: " << overlap << std::endl;
}

// Test case 2 to evaluate the overlap integral between an s-type and a p-type functions (l A = 0,  lB = 1) both centered at the origin
void test_case_2() {
    std::vector<double> centerA = {0.0, 0.0, 0.0};
    std::vector<double> centerB = {0.0, 0.0, 0.0};
    double alpha = 1.0;
    double beta = 1.0;

    // Create Shell objects for s-shell A and p-shell B
    Shell shellA(centerA, "s", alpha);
    Shell shellB(centerB, "p", beta);

    // Calculate the overlap integral as a vector of size 3
    std::vector<double> overlap = calculateOverlapIntegral_sp(shellA, shellB);

    std::cout << "Overlap integral between an s-shell and a p-shell is: " << std::endl;
    std::cout << "S^AB_spx = " << overlap[0] << std::endl;
    std::cout << "S^AB_spy = " << overlap[1] << std::endl;
    std::cout << "S^AB_spz = " << overlap[2] << std::endl;
}

// Test case 3 to evaluate the overlap integral between two p-type functions (l A = 1,  lB = 1) both centered at the origin
void test_case_3() {
    std::vector<double> centerA = {0.0, 0.0, 0.0};
    std::vector<double> centerB = {0.0, 0.0, 0.0};
    double alpha = 1.0;
    double beta = 1.0;

    // Create Shell objects for p-shells A and B
    Shell shellA(centerA, "p", alpha);
    Shell shellB(centerB, "p", beta);

    // Calculate the overlap integral as a matrix of size [3, 3]
    arma::mat overlap = calculateOverlapIntegral_pp(shellA, shellB);

    std::cout << "Overlap integral between two p-shells is: " << std::endl;
    std::cout << overlap(0, 0) << "   " << overlap(0, 1) << "   " << overlap(0, 2) << std::endl;
    std::cout << overlap(1, 0) << "   " << overlap(1, 1) << "   " << overlap(1, 2) << std::endl;
    std::cout << overlap(2, 0) << "   " << overlap(2, 1) << "   " << overlap(2, 2) << std::endl;
}

// after test cases with both shells centered at the origin
// a new function to take in two shells with different centers
// and parameters
void overlapIntegral(){
    // take user input for shell A
    double xA, yA, zA;
    cout << "Please enter the X coordinate for shell A:" << endl;
    cin >> xA;
    cout << "Please enter the Y coordinate for shell A:" << endl;
    cin >> yA;
    cout << "Please enter the Z coordinate for shell A:" << endl;
    cin >> zA;
    
    std::vector<double> centerA = {xA, yA, zA};

    string typeA;
    cout << "Please enter the type of shell A (s or p): " << endl;
    cin >> typeA;

    double alpha;
    cout << "Please enter the exponent of shell A: " << endl;
    cin >> alpha;

    // take user input for shell B
    double xB, yB, zB;
    cout << "Please enter the X coordinate for shell B:" << endl;
    cin >> xB;
    cout << "Please enter the Y coordinate for shell B:" << endl;
    cin >> yB;
    cout << "Please enter the Z coordinate for shell B:" << endl;
    cin >> zB;

    std::vector<double> centerB = {xB, yB, zB};

    string typeB;
    cout << "Please enter the type of shell B (s or p): " << endl;
    cin >> typeB;

    double beta;
    cout << "Please enter the exponent of shell B: " << endl;
    cin >> beta;

    // Create Shell objects for shell A and B
    // if typeA & typeB are both s, then shellA = typeA, shellB = typeB, calculateOverlapIntegral_ss
    // if typeA & typeB are both p, then shellA = typeA, shellB = typeB, calculateOverlapIntegral_pp
    // if typeA = s, typeB = p, then shellA = typeA, shellB = typeB, calculateOverlapIntegral_sp
    // if typeA = p, typeB = s, then shellA = typeB, shellB = typeA, calculateOverlapIntegral_sp
    if (typeA == "s" && typeB == "s") {
        Shell shellA(centerA, typeA, alpha);
        Shell shellB(centerB, typeB, beta);
        double overlap = calculateOverlapIntegral_ss(shellA, shellB);
        std::cout << "Overlap integral between two s-shells is: " << overlap << std::endl;
    } else if (typeA == "p" && typeB == "p") {
        Shell shellA(centerA, typeA, alpha);
        Shell shellB(centerB, typeB, beta);
        arma::mat overlap = calculateOverlapIntegral_pp(shellA, shellB);
        std::cout << "Overlap integral between two p-shells is: " << std::endl;
        std::cout << overlap(0, 0) << "   " << overlap(0, 1) << "   " << overlap(0, 2) << std::endl;
        std::cout << overlap(1, 0) << "   " << overlap(1, 1) << "   " << overlap(1, 2) << std::endl;
        std::cout << overlap(2, 0) << "   " << overlap(2, 1) << "   " << overlap(2, 2) << std::endl;
    } else if (typeA == "s" && typeB == "p") {
        Shell shellA(centerA, typeA, alpha);
        Shell shellB(centerB, typeB, beta);
        std::vector<double> overlap = calculateOverlapIntegral_sp(shellA, shellB);
        std::cout << "Overlap integral between an s-shell and a p-shell is: " << std::endl;
        std::cout << "S^AB_spx = " << overlap[0] << std::endl;
        std::cout << "S^AB_spy = " << overlap[1] << std::endl;
        std::cout << "S^AB_spz = " << overlap[2] << std::endl;
    } else if (typeA == "p" && typeB == "s") {
        Shell shellA(centerB, typeB, beta);
        Shell shellB(centerA, typeA, alpha);
        std::vector<double> overlap = calculateOverlapIntegral_sp(shellA, shellB);
        std::cout << "Overlap integral between an s-shell and a p-shell is: " << std::endl;
        std::cout << "S^AB_spx = " << overlap[0] << std::endl;
        std::cout << "S^AB_spy = " << overlap[1] << std::endl;
        std::cout << "S^AB_spz = " << overlap[2] << std::endl;
    }

}

int main() {

    test_case_1();
    test_case_2();
    test_case_3();
    overlapIntegral();


    return 0;
}
