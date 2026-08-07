#pragma once

#include <cstddef>

#if __has_include(<Eigen/Sparse>)
#include <Eigen/Sparse>
#elif __has_include(<eigen3/Eigen/Sparse>)
#include <eigen3/Eigen/Sparse>
#endif

Eigen::SparseMatrix<double> getSparseMassMatrix(
    std::size_t dim,
    const double uniteMass,
    const double radius,
    const int massNum
);
