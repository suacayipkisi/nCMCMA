#include "stdMatrixOperations.h"

//==========MULTIPLICATION================

std::vector<double> multbyValue_Vxval(const std::vector<double>& vec, const double value){
    const std::size_t T = vec.size();
    std::vector<double> resultVector = vec;
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        resultVector[i] *= value;
    }
    return resultVector;
}

std::vector<std::vector<double>> multbyValue_Mxval(const std::vector<std::vector<double>>& matrix, const double value){
    const std::size_t rows = matrix.size();
    if (rows == 0) return {};
    const std::size_t cols = matrix[0].size();
    std::vector<std::vector<double>> resultMatrix = matrix;
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < rows; ++i){
        for (std::size_t j = 0; j < cols; ++j){
            resultMatrix[i][j] *= value;
        }
    }
    return resultMatrix;
}

std::vector<std::vector<double>> mult_MxM(
    const std::vector<std::vector<double>>& matrix1,
    const std::vector<std::vector<double>>& matrix2
){
    const std::size_t R1 = matrix1.size();
    if (R1 == 0) return {};
    const std::size_t C1 = matrix1[0].size();
    const std::size_t C2 = matrix2[0].size();

    std::vector<std::vector<double>> resultMatrix(R1, std::vector<double>(C2, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R1; ++i){
        for (std::size_t j = 0; j < C2; ++j){
            for (std::size_t k = 0; k < C1; ++k){
                resultMatrix[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return resultMatrix;
}

std::vector<std::vector<double>> mult_VxVt(
    const std::vector<double>& vector1,
    const std::vector<double>& vector2
){
    const std::size_t T = vector1.size();
    const std::size_t U = vector2.size();
    std::vector<std::vector<double>> resultMatrix(T, std::vector<double>(U, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < U; ++j){
            resultMatrix[i][j] = vector1[i] * vector2[j];
        }
    }
    return resultMatrix;
}

std::vector<double> mult_MxV(
    const std::vector<std::vector<double>>& matrix1,
    const std::vector<double>& vector2
){
    const std::size_t R = matrix1.size();
    if (R == 0) return {};
    const std::size_t C = matrix1[0].size();
    std::vector<double> resultVector(R, 0.0);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R; ++i){
        for(std::size_t j = 0; j < C; ++j){
            resultVector[i] += matrix1[i][j] * vector2[j];
        }
    }
    return resultVector;
}

std::vector<std::vector<std::complex<double>>> multCompl_cMxcM(
    const std::vector<std::vector<std::complex<double>>>& complexMatrix1,
    const std::vector<std::vector<std::complex<double>>>& complexMatrix2
){
    const std::size_t R1 = complexMatrix1.size();
    if (R1 == 0) return {};
    const std::size_t C1 = complexMatrix1[0].size();
    const std::size_t C2 = complexMatrix2[0].size();

    std::vector<std::vector<std::complex<double>>> resultMatrix(R1, std::vector<std::complex<double>>(C2, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R1; ++i){
        for (std::size_t j = 0; j < C2; ++j){
            for (std::size_t k = 0; k < C1; ++k){
                resultMatrix[i][j] += complexMatrix1[i][k] * complexMatrix2[k][j];
            }
        }
    }
    return resultMatrix;
}

std::vector<std::vector<std::complex<double>>> multCompl_cMxM(
    const std::vector<std::vector<std::complex<double>>>& complexMatrix1,
    const std::vector<std::vector<double>>& matrix2
){
    const std::size_t R1 = complexMatrix1.size();
    if (R1 == 0) return {};
    const std::size_t C1 = complexMatrix1[0].size();
    const std::size_t C2 = matrix2[0].size();

    std::vector<std::vector<std::complex<double>>> resultMatrix(R1, std::vector<std::complex<double>>(C2, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R1; ++i){
        for (std::size_t j = 0; j < C2; ++j){
            for (std::size_t k = 0; k < C1; ++k){
                resultMatrix[i][j] += complexMatrix1[i][k] * matrix2[k][j];
            }
        }
    }
    return resultMatrix;
}

std::vector<std::vector<std::complex<double>>> multCompl_MxcM(
    const std::vector<std::vector<double>>& matrix1,
    const std::vector<std::vector<std::complex<double>>>& complexMatrix2
){
    const std::size_t R1 = matrix1.size();
    if (R1 == 0) return {};
    const std::size_t C1 = matrix1[0].size();
    const std::size_t C2 = complexMatrix2[0].size();

    std::vector<std::vector<std::complex<double>>> resultMatrix(R1, std::vector<std::complex<double>>(C2, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R1; ++i){
        for (std::size_t j = 0; j < C2; ++j){
            for (std::size_t k = 0; k < C1; ++k){
                resultMatrix[i][j] += matrix1[i][k] * complexMatrix2[k][j];
            }
        }
    }
    return resultMatrix;
}

std::vector<std::vector<std::complex<double>>> multCompl_MxcVal(
    const std::vector<std::vector<double>>& matrix1,
    const std::complex<double> value
){
    const std::size_t R = matrix1.size();
    if (R == 0) return {};
    const std::size_t C = matrix1[0].size();

    std::vector<std::vector<std::complex<double>>> resultMatrix(R, std::vector<std::complex<double>>(C, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R; ++i){
        for (std::size_t j = 0; j < C; ++j){
            resultMatrix[i][j] = matrix1[i][j] * value;
        }
    }
    return resultMatrix;
}

std::vector<std::vector<std::complex<double>>> multCompl_cMxVal(
    const std::vector<std::vector<std::complex<double>>>& matrix1,
    const double value
){
    const std::size_t R = matrix1.size();
    if (R == 0) return {};
    const std::size_t C = matrix1[0].size();

    std::vector<std::vector<std::complex<double>>> resultMatrix(R, std::vector<std::complex<double>>(C, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R; ++i){
        for (std::size_t j = 0; j < C; ++j){
            resultMatrix[i][j] = matrix1[i][j] * value;
        }
    }
    return resultMatrix;
}

std::vector<std::vector<std::complex<double>>> multCompl_cMxcVal(
    const std::vector<std::vector<std::complex<double>>>& matrix1,
    const std::complex<double> value
){
    const std::size_t R = matrix1.size();
    if (R == 0) return {};
    const std::size_t C = matrix1[0].size();

    std::vector<std::vector<std::complex<double>>> resultMatrix(R, std::vector<std::complex<double>>(C, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R; ++i){
        for (std::size_t j = 0; j < C; ++j){
            resultMatrix[i][j] = matrix1[i][j] * value;
        }
    }
    return resultMatrix;
}

//=============SUM===============

std::vector<double> sum_VxV(const std::vector<double>& vector1, const std::vector<double>& vector2){
    const std::size_t T = vector1.size();
    std::vector<double> resultMatrix(T, 0.0);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        resultMatrix[i] = vector1[i] + vector2[i];
    }
    return resultMatrix;
}

std::vector<std::vector<double>> sum_MxM(
    const std::vector<std::vector<double>>& matrix1,
    const std::vector<std::vector<double>>& matrix2
){
    const std::size_t R = matrix1.size();
    if (R == 0) return {};
    const std::size_t C = matrix1[0].size();

    std::vector<std::vector<double>> resultMatrix(R, std::vector<double>(C, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R; ++i){
        for (std::size_t j = 0; j < C; ++j){
            resultMatrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return resultMatrix;
}

std::vector<std::vector<std::complex<double>>> sumComplex_MxM(
    const std::vector<std::vector<std::complex<double>>>& matrix1,
    const std::vector<std::vector<std::complex<double>>>& matrix2
){
    const std::size_t R = matrix1.size();
    if (R == 0) return {};
    const std::size_t C = matrix1[0].size();

    std::vector<std::vector<std::complex<double>>> resultMatrix(R, std::vector<std::complex<double>>(C, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R; ++i){
        for (std::size_t j = 0; j < C; ++j){
            resultMatrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return resultMatrix;
}

std::vector<std::vector<std::complex<double>>> sumComplex_MxM(
    const std::vector<std::vector<double>>& matrix1,
    const std::vector<std::vector<std::complex<double>>>& matrix2
){
    const std::size_t R = matrix1.size();
    if (R == 0) return {};
    const std::size_t C = matrix1[0].size();

    std::vector<std::vector<std::complex<double>>> resultMatrix(R, std::vector<std::complex<double>>(C, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R; ++i){
        for (std::size_t j = 0; j < C; ++j){
            resultMatrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return resultMatrix;
}

std::vector<std::vector<std::complex<double>>> sumComplex_MxM(
    const std::vector<std::vector<std::complex<double>>>& matrix1,
    const std::vector<std::vector<double>>& matrix2
){
    const std::size_t R = matrix1.size();
    if (R == 0) return {};
    const std::size_t C = matrix1[0].size();

    std::vector<std::vector<std::complex<double>>> resultMatrix(R, std::vector<std::complex<double>>(C, 0.0));
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < R; ++i){
        for (std::size_t j = 0; j < C; ++j){
            resultMatrix[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }
    return resultMatrix;
}