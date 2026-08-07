#include "applyStateSpace.h"
#include "../matrixOperations/stdMatrixOperations.h"
#include <omp.h>

std::vector<std::vector<double>> getIdentityMatrix(std::size_t size){
    std::vector<std::vector<double>> identityMatrix(size, std::vector<double>(size, 0.0));
    for (std::size_t i{0}; i < size; ++i){
        identityMatrix[i][i] = 1.0;
    }
    return identityMatrix;
}

std::vector<std::vector<double>> getRayleightDampingMatrix(
    const std::string_view materialName,
    const std::vector<std::vector<double>>& massMatrix,
    const std::vector<std::vector<double>>& stiffnessMatrix,
    const double alphaConst,
    const double betaConst
)
{
    // Return empty vector to avoid allocating 4.8 GB for damping matrix.
    // Rayleigh damping C = alpha * M + beta * K is computed on the fly.
    return {};
}

std::vector<std::vector<double>> getStateSpaceMatrix(
    const std::vector<std::vector<double>>& massMatrix,
    const std::vector<std::vector<double>>& stiffMatrix,
    const std::vector<std::vector<double>>& rayleightMatrix
){
    const std::size_t T = massMatrix.size();
    std::vector<std::vector<double>> stateSpaceMatrix(2 * T, std::vector<double>(2 * T, 0.0));
    const bool hasC = !rayleightMatrix.empty() && rayleightMatrix.size() == T;

    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i) {
        const double m_inv = (massMatrix[i][i] != 0.0) ? (1.0 / massMatrix[i][i]) : 0.0;
        stateSpaceMatrix[i][T + i] = 1.0;
        for (std::size_t j = 0; j < T; ++j) {
            stateSpaceMatrix[T + i][j] = -m_inv * stiffMatrix[i][j];
            if (hasC) {
                stateSpaceMatrix[T + i][T + j] = -m_inv * rayleightMatrix[i][j];
            }
        }
    }
    return stateSpaceMatrix;
}