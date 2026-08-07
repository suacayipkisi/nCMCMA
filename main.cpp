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

#include "src/modalAnalysis/stiffMatrix.h"
#include "src/modalAnalysis/massMatrix.h"
#include "src/modalAnalysis/applyStateSpace.h"
#include "src/matrixOperations/stdEigenValueSolver.h"
#include "src/modalAnalysis/receptanceMatrix.h"

#ifndef MAIN_DIR
#define MAIN_DIR "main"
#endif

int main() {

    const unsigned int hardwareThreads{std::thread::hardware_concurrency()};
    const int numThreads{std::max(1, static_cast<int>((hardwareThreads > 0 ? hardwareThreads : 2) / 2))};
    omp_set_num_threads(numThreads);
    Eigen::setNbThreads(numThreads);

    std::cout << "Thread number in system: " << (hardwareThreads > 0 ? hardwareThreads : 1)
              << " | Using thread number: " << numThreads << "\n\n";

    constexpr std::size_t massNum{3};
    constexpr std::size_t dim{6 * massNum * massNum * massNum};
    
    std::string_view name{"Steel Spring A288"};

    const double alpha{.02287};
    const double beta{0.00000578};
    const double radius{0.01};
    const double density{7800};
    const double mass{4.0 / 3.0 * std::numbers::pi * radius * radius * radius};
    const double stiffnessConst{1000};

    static std::array<std::array<double, dim>, dim> massMatrix{getMassMatrix<dim>(mass, radius, static_cast<int>(massNum))};
    static std::array<std::array<double, dim>, dim> stiffnessMatrix{getStriffnessMatrix<dim>(stiffnessConst, radius, static_cast<int>(massNum))};
    static std::array<std::array<double, dim>, dim> rayleightDampingMatrix{getRayleightDampingMatrix<dim>(name, massMatrix, stiffnessMatrix)};

    std::array<std::array<double, 2 * dim>, 2 * dim>& stateSpaceMatrix{
        *new std::array<std::array<double, 2 * dim>, 2 * dim>{getStateSpaceMatrix<dim>(massMatrix, stiffnessMatrix, rayleightDampingMatrix)}
    };

    auto analysisResult{getAnalysisResult(stateSpaceMatrix)};

    const auto receptanceFrequency{std::real(analysisResult[2][1][1])};
    
    std::array<std::array<std::complex<double>, dim>, dim>& receptanceMatrix{
        *new std::array<std::array<std::complex<double>, dim>, dim>{getReceptanceMatrix(receptanceFrequency, massMatrix, stiffnessMatrix, rayleightDampingMatrix)}
    };

    std::string receptanceLocation{std::string(MAIN_DIR) + "/receptance_matrix.txt"};
    std::ofstream outFile(receptanceLocation);
    if (!outFile) {
        std::cerr << "[ERROR]: File cannot created!\n";
        delete &stateSpaceMatrix;
        delete &receptanceMatrix;
        return 1;
    }

    std::vector<std::string> rowBuffer(dim);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < dim; ++i) {
        std::ostringstream ss;
        for (std::size_t j = 0; j < dim; ++j) {
            ss << receptanceMatrix[i][j] << '\t';
        }
        ss << '\n';
        rowBuffer[i] = ss.str();
    }

    for (const auto& row : rowBuffer) {
        outFile << row;
    }

    std::cout << "[SUCCESS]Receptance matrix writed in 'receptance_matrix.txt' file succesfully.\n";

    std::array<std::complex<double>, dim> Q{};
    Q[0] = std::complex<double>(1000.0, 0.0);

    std::array<std::complex<double>, dim> q{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < dim; ++i) {
        std::complex<double> sum{0.0, 0.0};
        for (std::size_t j = 0; j < dim; ++j) {
            sum += receptanceMatrix[i][j] * Q[j];
        }
        q[i] = sum;
    }

    std::string displacementLocation{std::string(MAIN_DIR) + "/displacement.txt"};
    std::ofstream dispFile(displacementLocation);
    if (!dispFile) {
        std::cerr << "[ERROR]: displacement.txt file cannot be created!\n";
        delete &stateSpaceMatrix;
        delete &receptanceMatrix;
        return 1;
    }

    for (std::size_t i = 0; i < dim; ++i) {
        dispFile << "DOF " << i << ": " << q[i] << " (Magnitude: " << std::abs(q[i]) << ")\n";
    }

    std::cout << "[SUCCESS] Displacement vector written to 'displacement.txt' successfully.\n";

    delete &stateSpaceMatrix;
    delete &receptanceMatrix;

    return 0;
}