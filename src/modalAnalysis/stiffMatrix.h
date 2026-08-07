#pragma once

#include <vector>
#include <array>
#include <cstddef>

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

std::vector<std::vector<double>> getStriffnessMatrix(
    std::size_t dim,
    const double stiffConst,
    const double radius,
    const int massNum
);