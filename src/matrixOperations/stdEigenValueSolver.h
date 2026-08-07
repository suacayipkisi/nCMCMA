#pragma once

#if __has_include(<Eigen/Sparse>)
#include <Eigen/Sparse>
#elif __has_include(<eigen3/Eigen/Sparse>)
#include <eigen3/Eigen/Sparse>
#endif

double getNaturalFrequency(
    const Eigen::SparseMatrix<double>& M,
    const Eigen::SparseMatrix<double>& K,
    std::size_t modeIndex = 1
);
