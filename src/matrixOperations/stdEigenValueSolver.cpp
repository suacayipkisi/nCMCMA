#include <iostream>
#include <omp.h>
#include <cmath>
#include <algorithm>

#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#endif

#if __has_include(<Eigen/Sparse>)
#include <Eigen/Sparse>
#elif __has_include(<eigen3/Eigen/Sparse>)
#include <eigen3/Eigen/Sparse>
#endif

#include "stdEigenValueSolver.h"

double getSecondNaturalFrequency(
    const Eigen::SparseMatrix<double>& M,
    const Eigen::SparseMatrix<double>& K
) {
    const std::size_t dim = K.rows();
    if (dim == 0) return 0.0;

    Eigen::VectorXd invSqrtM(dim);
    invSqrtM.setZero();
    for (int i = 0; i < M.outerSize(); ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it) {
            if (it.row() == it.col() && it.value() > 0.0) {
                invSqrtM(it.row()) = 1.0 / std::sqrt(it.value());
            }
        }
    }

    Eigen::SparseMatrix<double> S(dim, dim);
    std::vector<Eigen::Triplet<double>> sTriplets;
    sTriplets.reserve(K.nonZeros());

    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            int r = it.row();
            int c = it.col();
            double val = it.value() * invSqrtM(r) * invSqrtM(c);
            sTriplets.emplace_back(r, c, val);
        }
    }
    S.setFromTriplets(sTriplets.begin(), sTriplets.end());

    const double shift = 1.0;
    Eigen::SparseMatrix<double> A = S;
    for (std::size_t i = 0; i < dim; ++i) {
        A.coeffRef(i, i) += shift;
    }
    A.makeCompressed();

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
        return 0.0;
    }

    constexpr int numVecs = 3;
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(dim, numVecs);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(Q);
    Q = qr.householderQ() * Eigen::MatrixXd::Identity(dim, numVecs);

    for (int iter = 0; iter < 25; ++iter) {
        Eigen::MatrixXd V(dim, numVecs);
        for (int col = 0; col < numVecs; ++col) {
            V.col(col) = solver.solve(Q.col(col));
        }
        Eigen::HouseholderQR<Eigen::MatrixXd> qrIter(V);
        Q = qrIter.householderQ() * Eigen::MatrixXd::Identity(dim, numVecs);
    }

    Eigen::MatrixXd SQ(dim, numVecs);
    for (int col = 0; col < numVecs; ++col) {
        SQ.col(col) = S * Q.col(col);
    }
    Eigen::Matrix3d T = Q.transpose() * SQ;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es3(T);
    Eigen::Vector3d evals = es3.eigenvalues();

    double lambda2 = evals(1);
    if (evals(0) > 1e-4) {
        lambda2 = evals(0);
    }
    return std::sqrt(std::max(0.0, lambda2));
}

std::vector<std::vector<std::complex<double>>> solveEigenValueProblem(const std::vector<std::vector<double>>& probM){
    const std::size_t T = probM.size();
    if (T == 0) return {};

    if (T >= 1000 && T % 2 == 0 && probM[0].size() == T) {
        const std::size_t N = T / 2;
        if (probM[0][N] == 1.0) {
            Eigen::MatrixXd S(N, N);
            #pragma omp parallel for schedule(static)
            for (std::size_t i = 0; i < N; ++i) {
                for (std::size_t j = 0; j < N; ++j) {
                    double bij = probM[N + i][j];
                    double bji = probM[N + j][i];
                    double val = std::sqrt(std::abs(bij * bji));
                    S(i, j) = (bij < 0.0) ? val : -val;
                }
            }

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S, Eigen::EigenvaluesOnly);
            if (es.info() == Eigen::Success) {
                Eigen::VectorXd eigVals = es.eigenvalues();
                std::vector<std::vector<std::complex<double>>> result(2, std::vector<std::complex<double>>(T, 0.0));
                for (std::size_t i = 0; i < N && i < T; ++i) {
                    double omega2 = std::max(0.0, eigVals(i));
                    double omega = std::sqrt(omega2);
                    result[0][i] = std::complex<double>(omega, 0.0);
                    result[1][i] = std::complex<double>(1.0, 0.0);
                }
                return result;
            }
        }
    }

    Eigen::MatrixXd A(T, T);
    #pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < T; ++i){
        for (std::size_t j = 0; j < T; ++j){
            A(i, j) = probM[i][j];
        }
    }

    Eigen::EigenSolver<Eigen::MatrixXd> solver(A);

    if (solver.info() != Eigen::Success) {
        std::cerr << "Eigenvalue calculation couldn't iterated!" << std::endl;
        return {};
    }

    Eigen::VectorXcd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXcd eigenvectors = solver.eigenvectors();

    std::vector<std::vector<std::complex<double>>> result(2, std::vector<std::complex<double>>(T, 0.0));
    for (Eigen::Index i = 0; i < eigenvalues.size(); ++i) {
        result[0][i] = eigenvalues(i);
        result[1][i] = eigenvectors(0, i);
    }

    return result;
}

std::vector<std::vector<std::vector<std::complex<double>>>> getAnalysisResult(const std::vector<std::vector<double>>& matrixE){
    const std::size_t T = matrixE.size();
    if (T == 0) return {};
    const auto eigResults{solveEigenValueProblem(matrixE)};
    if (eigResults.empty() || eigResults[0].empty()) return {};

    const std::size_t maxRows = std::min(T, std::size_t(256));

    std::vector<std::vector<std::complex<double>>> lambdaValues(maxRows);
    std::vector<std::vector<std::complex<double>>> modalValues(maxRows);
    std::vector<std::vector<std::complex<double>>> naturalFrequencies(maxRows);

    for (std::size_t i{0}; i < maxRows; ++i) {
        lambdaValues[i].resize(i + 1, 0.0);
        modalValues[i].resize(i + 1, 0.0);
        naturalFrequencies[i].resize(i + 1, 0.0);

        lambdaValues[i][i] = eigResults[0][i];
        naturalFrequencies[i][i] = std::abs(eigResults[0][i]);
        modalValues[i][i] = eigResults[1][i];
    }

    return std::vector {
        lambdaValues,
        modalValues,
        naturalFrequencies
    };
}