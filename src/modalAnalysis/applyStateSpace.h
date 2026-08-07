#pragma once

#include <array>
#include <cstddef>
#include <string_view>

#include "../matrixOperations/stdMatrixOperations.h"

template <std::size_t N>
std::array<std::array<double, N>, N> getIdentityMatrix(const size_t size){

    std::array<std::array<double, N>, N> identityMatrix;
    for (size_t i{0}; i< size; ++i){
        identityMatrix[i][i] = 1.0;
    }
    return identityMatrix;
}

template <std::size_t M>
std::array<std::array<double, M>, M> getRayleightDampingMatrix(
    const std::string_view materialName,
    const std::array<std::array<double, M>, M>& massMatrix,
    const std::array<std::array<double, M>, M>& stiffnessMatrix
)
{
    //normally we will pull these constants from file
    //Steel_Spring_A228
    double alphaConst {0.2287};
    double betaConst {0.00000578};

    static std::array<std::array<double, M>, M> rayleightDampingMatrix;
    rayleightDampingMatrix = sum_MxM(
        multbyValue_Mxval(massMatrix, alphaConst),
        multbyValue_Mxval(stiffnessMatrix, betaConst)
    );
    return rayleightDampingMatrix;
}

template <std::size_t T>
std::array<std::array<double, 2*T>, 2*T> getStateSpaceMatrix(
    const std::array<std::array<double, T>, T>& massMatrix,
    const std::array<std::array<double, T>, T>& stiffMatrix,
    const std::array<std::array<double, T>, T>& rayleightMatrix
){
    static std::array<std::array<double, T>, T> massInverseMatrix{massMatrix};
    for (size_t i{0}; i<T; ++i){
        if(massInverseMatrix[i][i] != 0)
        massInverseMatrix[i][i] = (1.0 / massInverseMatrix[i][i]);
    }
    static std::array<std::array<double, T>, T> min_Minv_K{
        multbyValue_Mxval(mult_MxM(massInverseMatrix, stiffMatrix),
        (-1.0))
    };
    static std::array<std::array<double, T>, T> min_Minv_C{
        multbyValue_Mxval(mult_MxM(massInverseMatrix, rayleightMatrix),
        (-1.0))
    };

    static std::array<std::array<double, 2*T>, 2*T> stateSpaceMatrix{};
    for (size_t i{0}; i<T; ++i){
        for (size_t j{0}; j<T; ++j){
            if(i==j){
                stateSpaceMatrix[i][T+j] = 1.0;
            }
            stateSpaceMatrix[T+i][j] = min_Minv_K[i][j];
            stateSpaceMatrix[T+i][T+j] = min_Minv_C[i][j];
        }
    }
    return stateSpaceMatrix;
}