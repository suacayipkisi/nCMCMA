#pragma once

#include <complex>
#include <vector>

#if __has_include(<Eigen/Sparse>)
#include <Eigen/Sparse>
#elif __has_include(<eigen3/Eigen/Sparse>)
#include <eigen3/Eigen/Sparse>
#endif

std::vector<std::vector<std::complex<double>>> solveEigenValueProblem(const std::vector<std::vector<double>>& probM);

std::vector<std::vector<std::vector<std::complex<double>>>> getAnalysisResult(const std::vector<std::vector<double>>& matrixE);

double getNaturalFrequency(
    const Eigen::SparseMatrix<double>& M,
    const Eigen::SparseMatrix<double>& K,
    std::size_t modeIndex = 1
);
