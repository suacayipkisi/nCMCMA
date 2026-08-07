#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <array>
#include <complex>
#include <cstddef>
#include <thread>
#include <algorithm>
#include <omp.h>
#include <numbers>

#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#else
#error "Eigen headers not found. Add Eigen include path (e.g., /usr/include/eigen3)."
#endif

#include "src/modalAnalysis/stiffMatrix.h"
#include "src/modalAnalysis/massMatrix.h"
#include "src/matrixOperations/stdEigenValueSolver.h"

#ifndef MAIN_DIR
#define MAIN_DIR "main"
#endif

int main() {

    const unsigned int hardwareThreads{std::thread::hardware_concurrency()};
    const int numThreads{static_cast<int>(hardwareThreads > 0 ? hardwareThreads : 4)};
    omp_set_num_threads(numThreads);
    Eigen::setNbThreads(1); // Eigen inner parallelism off to avoid OpenMP thread collision

    std::cout << "Thread number in system: " << (hardwareThreads > 0 ? hardwareThreads : 1)
              << " | Using thread number: " << numThreads << "\n\n";

    // defining the cube size
    constexpr std::size_t massNum{16};
    constexpr std::size_t dim{6 * massNum * massNum * massNum};
    
    std::string_view name{"Steel Spring A288"};

    // defining material properties
    //normally we will pull these constants from file
    //Steel_Spring_A228
    const double alpha{.02287};
    const double beta{0.00000578};
    const double radius{0.01};
    const double density{7800};
    const double mass{4.0 / 3.0 * std::numbers::pi * radius * radius * radius};
    const double stiffnessConst{1000};

    // Build sparse mass and stiffness matrices (uses ~12 MB RAM instead of 9.6 GB)
    std::cout << "Building sparse mass and stiffness matrices...\n";
    Eigen::SparseMatrix<double> sparseM = getSparseMassMatrix(dim, mass, radius, static_cast<int>(massNum));
    Eigen::SparseMatrix<double> sparseK = getSparseStiffnessMatrix(dim, stiffnessConst, radius, static_cast<int>(massNum));

    std::cout << "Calculating natural frequency using sparse eigenvalue solver...\n";

    const double receptanceFrequency = getNaturalFrequency(sparseM, sparseK, 2);
    std::cout << "Calculated natural frequency: " << receptanceFrequency << " rad/s\n\n";

    const double w = receptanceFrequency;
    const double w2 = w * w;

    // Build sparse complex dynamic stiffness matrix Z = K - w^2 * M + j * w * (alpha * M + beta * K)
    std::complex<double> k_coeff(1.0, w * beta);
    std::complex<double> m_coeff(-w2, w * alpha);

    Eigen::SparseMatrix<std::complex<double>> Z(dim, dim);
    std::vector<Eigen::Triplet<std::complex<double>>> zTriplets;
    zTriplets.reserve(sparseK.nonZeros() + sparseM.nonZeros());

    for (int k = 0; k < sparseK.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(sparseK, k); it; ++it) {
            zTriplets.emplace_back(it.row(), it.col(), std::complex<double>(it.value(), 0.0) * k_coeff);
        }
    }
    for (int k = 0; k < sparseM.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(sparseM, k); it; ++it) {
            zTriplets.emplace_back(it.row(), it.col(), std::complex<double>(it.value(), 0.0) * m_coeff);
        }
    }
    Z.setFromTriplets(zTriplets.begin(), zTriplets.end());

    std::cout << "Factoring sparse impedance matrix Z (" << dim << "x" << dim << ")...\n";
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> lu;
    lu.compute(Z);
    if (lu.info() != Eigen::Success) {
        std::cerr << "[ERROR]: SparseLU factorization failed!\n";
        return 1;
    }

    // 1. Solve displacement vector q = Z^{-1} * Q
    std::cout << "Solving displacement vector q...\n";
    Eigen::VectorXcd Q_vec = Eigen::VectorXcd::Zero(dim);
    Q_vec(0) = std::complex<double>(1000.0, 0.0);
    Eigen::VectorXcd q_vec = lu.solve(Q_vec);

    std::string displacementLocation{std::string(MAIN_DIR) + "/displacement.txt"};
    std::ofstream dispFile(displacementLocation);
    if (!dispFile) {
        std::cerr << "[ERROR]: displacement.txt file cannot be created!\n";
        return 1;
    }

    for (std::size_t i = 0; i < dim; ++i) {
        dispFile << "DOF " << i << ": " << q_vec(i) << " (Magnitude: " << std::abs(q_vec(i)) << ")\n";
    }
    std::cout << "[SUCCESS] Displacement vector written to 'displacement.txt' successfully.\n";

    return 0;
}