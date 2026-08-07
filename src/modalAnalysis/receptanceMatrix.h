#pragma once

#include <array>
#include <complex>
#include <iostream>
#include <stdexcept>

#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#else
#error "Eigen headers not found. Add Eigen include path (e.g., /usr/include/eigen3)."
#endif

#include "../matrixOperations/stdMatrixOperations.h"

template <std::size_t T>
std::array<std::array<std::complex<double>, T>, T> getReceptanceMatrix(
    const double givenFrequency,
    const std::array<std::array<double, T>, T>& massMatrix,
    const std::array<std::array<double, T>, T>& stiffnessMatrix,
    const std::array<std::array<double, T>, T>& rayleightDampingMatrix
){
    static std::array<std::array<double, T>, T> min_w2_M;
    min_w2_M = multbyValue_Mxval(massMatrix, ((-1.0) * givenFrequency * givenFrequency));
    static std::array<std::array<std::complex<double>, T>, T> j_w_C;
    j_w_C = multCompl_MxcVal(rayleightDampingMatrix, (std::complex<double>(0, 1) * givenFrequency));
    static std::array<std::array<std::complex<double>, T>, T> rMatrix;
    rMatrix = sumComplex_MxM(sumComplex_MxM(min_w2_M, j_w_C), stiffnessMatrix);

    Eigen::MatrixXcd systemMatrix(T, T);
    systemMatrix.setZero();
    for (std::size_t i = 0; i < T; ++i) {
        for (std::size_t j = 0; j < T; ++j) {
            systemMatrix(i, j) = rMatrix[i][j];
        }
    }

    Eigen::MatrixXcd identity = Eigen::MatrixXcd::Identity(T, T);
    Eigen::FullPivLU<Eigen::MatrixXcd> lu(systemMatrix);
    if (lu.rank() < T) {
        throw std::runtime_error("[ERROR]: Matrix is singular and cannot be inverted.");
    }

    Eigen::MatrixXcd inverseMatrix = lu.solve(identity);

    std::array<std::array<std::complex<double>, T>, T> receptanceMatrix{};
    for (std::size_t i = 0; i < T; ++i) {
        for (std::size_t j = 0; j < T; ++j) {
            receptanceMatrix[i][j] = inverseMatrix(i, j);
        }
    }

    return receptanceMatrix;
}