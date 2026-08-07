#pragma once

#include <cstddef>
#include <array>
#include <complex>
#include <omp.h>

//==========MULTİPLİCATİON================

//value * vector or vector * value
template <std::size_t T>
std::array<double, T> multbyValue_Vxval(const std::array<double, T>& vector, const double value){
    std::array<double, T> resultVector{vector};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        resultVector[i] *= value;
    }
    return resultVector;
}

//matrix * value or value * matrix
template <std::size_t T>
std::array<std::array<double, T>, T> multbyValue_Mxval(const std::array<std::array<double, T>, T>& matrix, const double value){
    std::array<std::array<double, T>, T> resultMatrix{matrix};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < T; ++j){
            resultMatrix[i][j] *= value;
        }
    }
    return resultMatrix;
}

//matrix * matrix
template<std::size_t T, std::size_t U, std::size_t V>
std::array<std::array<double, T>, V> mult_MxM(const std::array<std::array<double, T>, U>& matrix1, const std::array<std::array<double, U>, V>& matrix2){
    std::array<std::array<double, T>, V> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < V; ++j){
            for (std::size_t k = 0; k < U; ++k){
                resultMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return resultMatrix;
}


//vector * vector_transpose
template<std::size_t T, std::size_t U>
std::array<std::array<double, T>, U> mult_VxVt(const std::array<double, T>& vector1, const std::array<double, U>& matrix2){
    std::array<std::array<double, T>, U> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < U; ++j){
            resultMatrix[i][j] = vector1[i] * matrix2[j];
        }
    }
    return resultMatrix;
}

//matrix * vector
template<std::size_t T, std::size_t U>
std::array<double, T> mult_MxV(const std::array<std::array<double, T>, U>& matrix1, const std::array<double, U>& vector2){
    std::array<double, T> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for(std::size_t j = 0; j < U; ++j){
            resultMatrix[i] += matrix1[i][j] * vector2[j];
        }
    }
    return resultMatrix;
}

// complexMatrix * complexMatrix
template<std::size_t T, std::size_t U, std::size_t V>
std::array<std::array<std::complex<double>, T>, V> multCompl_cMxcM(
    const std::array<std::array<std::complex<double>, T>, U>& complexMatrix1,
    const std::array<std::array<std::complex<double>, U>, V>& complexMatrix2
){
    std::array<std::array<std::complex<double>, T>, V> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < V; ++j){
            for (std::size_t k = 0; k < U; ++k){
                resultMatrix[i][j] += complexMatrix1[i][k] * complexMatrix2[k][j];
            }
        }
    }
    return resultMatrix;
}

// complexMatrix * Matrix
template<std::size_t T, std::size_t U, std::size_t V>
std::array<std::array<std::complex<double>, T>, V> multCompl_cMxM(
    const std::array<std::array<std::complex<double>, T>, U>& complexMatrix1,
    const std::array<std::array<double, U>, V>& matrix2
){
    std::array<std::array<std::complex<double>, T>, V> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < V; ++j){
            for (std::size_t k = 0; k < U; ++k){
                resultMatrix[i][j] += complexMatrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return resultMatrix;
}

// Matrix * complexMatrix
template<std::size_t T, std::size_t U, std::size_t V>
std::array<std::array<std::complex<double>, T>, V> multCompl_MxcM(
    const std::array<std::array<double, T>, U>& matrix1,
    const std::array<std::array<std::complex<double>, U>, V>& complexMatrix2
){
    std::array<std::array<std::complex<double>, T>, V> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < V; ++j){
            for (std::size_t k = 0; k < U; ++k){
                resultMatrix[i][j] += matrix1[i][k] * complexMatrix2[k][j];
            }
        }
    }
    return resultMatrix;
}

// matrix * complexVal
template<std::size_t T, std::size_t U>
std::array<std::array<std::complex<double>, T>, U> multCompl_MxcVal(
    const std::array<std::array<double, T>, U>& matrix1,
    const std::complex<double> value
){
    std::array<std::array<std::complex<double>, T>, U> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < U; ++j){
            resultMatrix[i][j] = matrix1[i][j] * value;
        }
    }
    return resultMatrix;
}

// complexMatrix * Val
template<std::size_t T, std::size_t U>
std::array<std::array<std::complex<double>, T>, U> multCompl_cMxVal(
    const std::array<std::array<std::complex<double>, T>, U>& matrix1,
    const double value
){
    std::array<std::array<std::complex<double>, T>, U> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < U; ++j){
            resultMatrix[i][j] = matrix1[i][j] * value;
        }
    }
    return resultMatrix;
}

// complexMatrix * complexValue
template<std::size_t T, std::size_t U>
std::array<std::array<std::complex<double>, T>, U> multCompl_cMxcVal(
    const std::array<std::array<std::complex<double>, T>, U>& matrix1,
    const std::complex<double> value
){
    std::array<std::array<std::complex<double>, T>, U> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < U; ++j){
            resultMatrix[i][j] = matrix1[i][j] * value;
        }
    }
    return resultMatrix;
}


//=============SUM===============

// vec + vec
template <std::size_t T>
std::array<double, T> sum_VxV(const std::array<double, T>& vector1, const std::array<double, T>& vector2){
    std::array<double, T> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        resultMatrix[i] = vector1[i] + vector2[i];
    }
    return resultMatrix;
}

// matrix + matrix
template <std::size_t T, std::size_t U>
std::array<std::array<double, T>, U> sum_MxM(
    const std::array<std::array<double, T>, U>& matrix1,
    const std::array<std::array<double, T>, U>& matrix2
){
    std::array<std::array<double, T>, U> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < U; ++j){
            resultMatrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return resultMatrix;
}

// complexMatrix + complexMatrix

template <std::size_t T, std::size_t U>
std::array<std::array<std::complex<double>, T>, U> sumComplex_MxM(
    const std::array<std::array<std::complex<double>, T>, U>& matrix1,
    const std::array<std::array<std::complex<double>, T>, U>& matrix2
){
    std::array<std::array<std::complex<double>, T>, U> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < U; ++j){
            resultMatrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return resultMatrix;
}
// matrix + complexMatrix

template <std::size_t T, std::size_t U>
std::array<std::array<std::complex<double>, T>, U> sumComplex_MxM(
    const std::array<std::array<double, T>, U>& matrix1,
    const std::array<std::array<std::complex<double>, T>, U>& matrix2
){
    std::array<std::array<std::complex<double>, T>, U> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < U; ++j){
            resultMatrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return resultMatrix;
}
//complexMatrix + matrix

template <std::size_t T, std::size_t U>
std::array<std::array<std::complex<double>, T>, U> sumComplex_MxM(
    const std::array<std::array<std::complex<double>, T>, U>& matrix1,
    const std::array<std::array<double, T>, U>& matrix2
){
    std::array<std::array<std::complex<double>, T>, U> resultMatrix{};
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < U; ++j){
            resultMatrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return resultMatrix;
}