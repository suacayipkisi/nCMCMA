#pragma once

#include <cstddef>
#include <vector>
#include <complex>

//==========MULTIPLICATION================

//value * vector or vector * value
std::vector<double> multbyValue_Vxval(const std::vector<double>& vec, const double value);

//matrix * value or value * matrix
std::vector<std::vector<double>> multbyValue_Mxval(const std::vector<std::vector<double>>& matrix, const double value);

//matrix * matrix
std::vector<std::vector<double>> mult_MxM(
    const std::vector<std::vector<double>>& matrix1,
    const std::vector<std::vector<double>>& matrix2
);

//vector * vector_transpose
std::vector<std::vector<double>> mult_VxVt(
    const std::vector<double>& vector1,
    const std::vector<double>& vector2
);

//matrix * vector
std::vector<double> mult_MxV(
    const std::vector<std::vector<double>>& matrix1,
    const std::vector<double>& vector2
);

// complexMatrix * complexMatrix
std::vector<std::vector<std::complex<double>>> multCompl_cMxcM(
    const std::vector<std::vector<std::complex<double>>>& complexMatrix1,
    const std::vector<std::vector<std::complex<double>>>& complexMatrix2
);

// complexMatrix * Matrix
std::vector<std::vector<std::complex<double>>> multCompl_cMxM(
    const std::vector<std::vector<std::complex<double>>>& complexMatrix1,
    const std::vector<std::vector<double>>& matrix2
);

// Matrix * complexMatrix
std::vector<std::vector<std::complex<double>>> multCompl_MxcM(
    const std::vector<std::vector<double>>& matrix1,
    const std::vector<std::vector<std::complex<double>>>& complexMatrix2
);

// matrix * complexVal
std::vector<std::vector<std::complex<double>>> multCompl_MxcVal(
    const std::vector<std::vector<double>>& matrix1,
    const std::complex<double> value
);

// complexMatrix * Val
std::vector<std::vector<std::complex<double>>> multCompl_cMxVal(
    const std::vector<std::vector<std::complex<double>>>& matrix1,
    const double value
);

// complexMatrix * complexValue
std::vector<std::vector<std::complex<double>>> multCompl_cMxcVal(
    const std::vector<std::vector<std::complex<double>>>& matrix1,
    const std::complex<double> value
);

//=============SUM===============

// vec + vec
std::vector<double> sum_VxV(const std::vector<double>& vector1, const std::vector<double>& vector2);

// matrix + matrix
std::vector<std::vector<double>> sum_MxM(
    const std::vector<std::vector<double>>& matrix1,
    const std::vector<std::vector<double>>& matrix2
);

// complexMatrix + complexMatrix
std::vector<std::vector<std::complex<double>>> sumComplex_MxM(
    const std::vector<std::vector<std::complex<double>>>& matrix1,
    const std::vector<std::vector<std::complex<double>>>& matrix2
);

// matrix + complexMatrix
std::vector<std::vector<std::complex<double>>> sumComplex_MxM(
    const std::vector<std::vector<double>>& matrix1,
    const std::vector<std::vector<std::complex<double>>>& matrix2
);

//complexMatrix + matrix
std::vector<std::vector<std::complex<double>>> sumComplex_MxM(
    const std::vector<std::vector<std::complex<double>>>& matrix1,
    const std::vector<std::vector<double>>& matrix2
);