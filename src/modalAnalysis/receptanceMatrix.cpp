#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#else
#error "Eigen headers not found. Add Eigen include path (e.g., /usr/include/eigen3)."
#endif

#include "receptanceMatrix.h"
#include "../matrixOperations/stdMatrixOperations.h"
#include <omp.h>
#include <stdexcept>

std::vector<std::vector<std::complex<double>>> getReceptanceMatrix(
    const double givenFrequency,
    const std::vector<std::vector<double>>& massMatrix,
    const std::vector<std::vector<double>>& stiffnessMatrix,
    const std::vector<std::vector<double>>& rayleightDampingMatrix
){
    const std::size_t T = massMatrix.size();
    const double w2 = givenFrequency * givenFrequency;
    const std::complex<double> j_w(0.0, givenFrequency);

    const double alpha = 0.02287;
    const double beta = 0.00000578;
    const bool hasC = !rayleightDampingMatrix.empty() && rayleightDampingMatrix.size() == T;

    Eigen::MatrixXcd systemMatrix(T, T);

    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i) {
        for (std::size_t j = 0; j < T; ++j) {
            double m_val = (i == j) ? massMatrix[i][i] : 0.0;
            double c_val = hasC ? rayleightDampingMatrix[i][j] : (m_val * alpha + stiffnessMatrix[i][j] * beta);
            systemMatrix(i, j) = stiffnessMatrix[i][j] - w2 * m_val + j_w * c_val;
        }
    }

    Eigen::PartialPivLU<Eigen::MatrixXcd> lu(systemMatrix);
    systemMatrix = Eigen::MatrixXcd(); // Free 9.6 GB RAM immediately!

    Eigen::MatrixXcd inverseMatrix = lu.solve(Eigen::MatrixXcd::Identity(T, T));

    std::vector<std::vector<std::complex<double>>> receptanceMatrix(T, std::vector<std::complex<double>>(T, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i) {
        for (std::size_t j = 0; j < T; ++j) {
            receptanceMatrix[i][j] = inverseMatrix(i, j);
        }
    }

    inverseMatrix = Eigen::MatrixXcd(); // Free 9.6 GB RAM immediately!

    return receptanceMatrix;
}