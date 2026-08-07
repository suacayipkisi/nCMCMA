#pragma once

#include <vector>
#include <complex>

std::vector<std::vector<std::complex<double>>> getReceptanceMatrix(
    const double givenFrequency,
    const std::vector<std::vector<double>>& massMatrix,
    const std::vector<std::vector<double>>& stiffnessMatrix,
    const std::vector<std::vector<double>>& rayleightDampingMatrix
);