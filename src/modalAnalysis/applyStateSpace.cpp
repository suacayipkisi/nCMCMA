#include "applyStateSpace.h"
#include "../matrixOperations/stdMatrixOperations.h"

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
    return sum_MxM(
        multbyValue_Mxval(massMatrix, alphaConst),
        multbyValue_Mxval(stiffnessMatrix, betaConst)
    );
}

std::vector<std::vector<double>> getStateSpaceMatrix(
    const std::vector<std::vector<double>>& massMatrix,
    const std::vector<std::vector<double>>& stiffMatrix,
    const std::vector<std::vector<double>>& rayleightMatrix
){
    const std::size_t T = massMatrix.size();
    std::vector<std::vector<double>> massInverseMatrix = massMatrix;
    for (std::size_t i{0}; i < T; ++i){
        if (massInverseMatrix[i][i] != 0)
            massInverseMatrix[i][i] = (1.0 / massInverseMatrix[i][i]);
    }
    std::vector<std::vector<double>> min_Minv_K = 
        multbyValue_Mxval(mult_MxM(massInverseMatrix, stiffMatrix), (-1.0));

    std::vector<std::vector<double>> min_Minv_C = 
        multbyValue_Mxval(mult_MxM(massInverseMatrix, rayleightMatrix), (-1.0));

    std::vector<std::vector<double>> stateSpaceMatrix(2 * T, std::vector<double>(2 * T, 0.0));
    for (std::size_t i{0}; i < T; ++i){
        for (std::size_t j{0}; j < T; ++j){
            if (i == j){
                stateSpaceMatrix[i][T + j] = 1.0;
            }
            stateSpaceMatrix[T + i][j] = min_Minv_K[i][j];
            stateSpaceMatrix[T + i][T + j] = min_Minv_C[i][j];
        }
    }
    return stateSpaceMatrix;
}