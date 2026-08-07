#pragma once

#include <vector>
#include <cstddef>

std::vector<std::vector<double>> getMassMatrix(
    std::size_t dim,
    const double uniteMass,
    const double radius,
    const int massNum
);
