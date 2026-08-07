#pragma once

#include <vector>
#include <cstddef>
#include <string_view>

std::vector<std::vector<double>> getIdentityMatrix(std::size_t size);

std::vector<std::vector<double>> getRayleightDampingMatrix(
    const std::string_view materialName,
    const std::vector<std::vector<double>>& massMatrix,
    const std::vector<std::vector<double>>& stiffnessMatrix,
    const double alphaConst,
    const double betaConst
);

std::vector<std::vector<double>> getStateSpaceMatrix(
    const std::vector<std::vector<double>>& massMatrix,
    const std::vector<std::vector<double>>& stiffMatrix,
    const std::vector<std::vector<double>>& rayleightMatrix
);