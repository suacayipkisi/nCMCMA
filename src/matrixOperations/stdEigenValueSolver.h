#pragma once

#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#else
#error "Eigen headers not found. Add Eigen include path (e.g., /usr/include/eigen3)."
#endif

#include <complex>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>

template <std::size_t T>
std::vector<std::vector<std::complex<double>>> solveEigenValueProblem(const std::array<std::array<double, T>, T>& probM){
    Eigen::MatrixXd A(T, T);
    for (std::size_t i{0}; i < T; ++i){
        for (std::size_t j{0}; j < T; ++j){
            A(i, j) = probM[i][j];
        }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> solver(A);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue calculation couldn't iterated!" << std::endl;
        return {};
    }

    Eigen::VectorXcd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXcd eigenvectors = solver.eigenvectors();

    std::vector<std::vector<std::complex<double>>> result(2, std::vector<std::complex<double>>(T, 0.0));
    for (Eigen::Index i = 0; i < eigenvalues.size(); ++i) {
        result[0][i] = eigenvalues(i);
        result[1][i] = eigenvectors(0, i);
    }

    return result;
}

template <std::size_t T>
std::vector<std::vector<std::vector<std::complex<double>>>> getAnalysisResult(const std::array<std::array<double, T>, T>& matrixE){
    const auto eigResults{solveEigenValueProblem(matrixE)};

    std::vector<std::vector<std::complex<double>>> lambdaValues(T, std::vector<std::complex<double>>(T, 0.0));
    std::vector<std::vector<std::complex<double>>> modalValues(T, std::vector<std::complex<double>>(T, 0.0));
    std::vector<std::vector<std::complex<double>>> naturalFrequencies(T, std::vector<std::complex<double>>(T, 0.0));

    for (std::size_t i{0}; i < T; ++i) {
        lambdaValues[i][i] = eigResults[0][i];
        naturalFrequencies[i][i] = std::abs(eigResults[0][i]);
        modalValues[i][i] = eigResults[1][i];
    }

    return std::vector {
        lambdaValues,
        modalValues,
        naturalFrequencies
    };
}