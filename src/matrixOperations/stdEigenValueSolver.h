#pragma once

#include <complex>
#include <vector>

std::vector<std::vector<std::complex<double>>> solveEigenValueProblem(const std::vector<std::vector<double>>& probM);

std::vector<std::vector<std::vector<std::complex<double>>>> getAnalysisResult(const std::vector<std::vector<double>>& matrixE);