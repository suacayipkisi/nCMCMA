#pragma once

#include <array>
#include <cstddef>

template <size_t T>
std::array<std::array<double, T>, T> getMassMatrix(
    const double uniteMass,
    const double radius,
    const int massNum //real mass num = massNum^3
){
    const int totalMassNum{massNum * massNum * massNum};
    static std::array<std::array<double, T>, T> massMatrix{};
    for (int i{0}; i < totalMassNum; ++i){
        const int base{6 * i};
        massMatrix[base + 0][base + 0] = uniteMass;
        massMatrix[base + 1][base + 1] = uniteMass;
        massMatrix[base + 2][base + 2] = uniteMass;
        massMatrix[base + 3][base + 3] = uniteMass;
        massMatrix[base + 4][base + 4] = uniteMass * (2.0 / 5.0) * radius * radius;
        massMatrix[base + 5][base + 5] = uniteMass * (2.0 / 5.0) * radius * radius;
    }
    return massMatrix;
}
