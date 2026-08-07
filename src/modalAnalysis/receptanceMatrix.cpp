#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#else
#error "Eigen headers not found. Add Eigen include path (e.g., /usr/include/eigen3)."
#endif

#include "receptanceMatrix.h"
#include "../matrixOperations/stdMatrixOperations.h"

std::vector<std::vector<std::complex<double>>> getReceptanceMatrix(
    const double givenFrequency,
    const std::vector<std::vector<double>>& massMatrix,
    const std::vector<std::vector<double>>& stiffnessMatrix,
    const std::vector<std::vector<double>>& rayleightDampingMatrix
){
    const std::size_t T = massMatrix.size();
    std::vector<std::vector<double>> min_w2_M = 
        multbyValue_Mxval(massMatrix, ((-1.0) * givenFrequency * givenFrequency));
    std::vector<std::vector<std::complex<double>>> j_w_C = 
        multCompl_MxcVal(rayleightDampingMatrix, (std::complex<double>(0, 1) * givenFrequency));
    std::vector<std::vector<std::complex<double>>> rMatrix = 
        sumComplex_MxM(sumComplex_MxM(min_w2_M, j_w_C), stiffnessMatrix);

    Eigen::MatrixXcd systemMatrix(T, T);
    systemMatrix.setZero();
    for (std::size_t i = 0; i < T; ++i) {
        for (std::size_t j = 0; j < T; ++j) {
            systemMatrix(i, j) = rMatrix[i][j];
        }
    }

    Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(T, T);
    Eigen::FullPivLU<Eigen::MatrixXcd> lu(systemMatrix);
    if (lu.rank() < static_cast<Eigen::Index>(T)) {
        throw std::runtime_error("[ERROR]: Matrix is singular and cannot be inverted.");
    }

    Eigen::MatrixXcd inverseMatrix = lu.solve(identity);

    std::vector<std::vector<std::complex<double>>> receptanceMatrix(T, std::vector<std::complex<double>>(T, 0.0));
    for (std::size_t i = 0; i < T; ++i) {
        for (std::size_t j = 0; j < T; ++j) {
            receptanceMatrix[i][j] = inverseMatrix(i, j);
        }
    }

    return receptanceMatrix;
}