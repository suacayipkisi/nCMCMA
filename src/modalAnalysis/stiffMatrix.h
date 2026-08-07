#pragma once

#include <vector>
#include <array>
#include <cstddef>

#if __has_include(<Eigen/Sparse>)
#include <Eigen/Sparse>
#elif __has_include(<eigen3/Eigen/Sparse>)
#include <eigen3/Eigen/Sparse>
#endif

void fillDiag(
    const std::array<int, 3>& elementPosition,
    const int elementID,
    std::vector<std::vector<double>>& stiffnessMatrix,
    const double stiffConst,
    const double radius,
    const int massNum
);

void fillElement(
    const std::array<int, 3>& elementPosition,
    const int elementID,
    const std::array<int, 3>& effectorPosition,
    const int effectorID,
    std::vector<std::vector<double>>& stiffnessMatrix,
    const double stiffConst,
    const double radius
);

std::vector<std::vector<double>> getStiffnessMatrix(
    std::size_t dim,
    const double stiffConst,
    const double radius,
    const int massNum
);

Eigen::SparseMatrix<double> getSparseStiffnessMatrix(
    std::size_t dim,
    const double stiffConst,
    const double radius,
    const int massNum
);