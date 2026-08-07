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

double getNaturalFrequency(
    const Eigen::SparseMatrix<double>& M,
    const Eigen::SparseMatrix<double>& K,
    std::size_t modeIndex
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

    const int numVecs = static_cast<int>(std::max(std::size_t(3), modeIndex + 2));
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(dim, numVecs);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(Q);
    Q = qr.householderQ() * Eigen::MatrixXd::Identity(dim, numVecs);

    for (int iter = 0; iter < 25; ++iter) {
        Eigen::MatrixXd V(dim, numVecs);
        #pragma omp parallel for schedule(static)
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
    Eigen::MatrixXd T = Q.transpose() * SQ;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
    Eigen::VectorXd evals = es.eigenvalues();

    std::size_t idx = std::min(modeIndex, static_cast<std::size_t>(evals.size() - 1));
    double lambda = evals(static_cast<Eigen::Index>(idx));

    return std::sqrt(std::max(0.0, lambda));
}
