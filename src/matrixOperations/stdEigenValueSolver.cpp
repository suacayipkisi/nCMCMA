#include <Eigen/Core>
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
)
{
    const std::size_t dim = K.rows();
    if (dim == 0 || M.rows() != dim) return 0.0;

    Eigen::VectorXd invSqrtM = Eigen::VectorXd::Zero(dim);
    for (int i = 0; i < M.outerSize(); ++i) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it) {
            if (it.row() == it.col() && it.value() > 0.0) {
                invSqrtM(it.row()) = 1.0 / std::sqrt(it.value());
            }
        }
    }

    std::vector<Eigen::Triplet<double>> sTriplets;
    sTriplets.reserve(K.nonZeros());

    for (int k = 0; k < K.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(K, k); it; ++it) {
            Eigen::Index r = it.row();
            Eigen::Index c = it.col();
            if(r < dim && c < dim) {
                double val = it.value() * invSqrtM(r) * invSqrtM(c);
                sTriplets.emplace_back(r, c, val);
            }    
        }
    }
    Eigen::SparseMatrix<double> S(dim, dim);
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
    
    if (static_cast<Eigen::Index>(modeIndex) >= dim) {
        return 0.0; 
    }

    const Eigen::Index requestedVecs = static_cast<Eigen::Index>(modeIndex) + 4;
    const Eigen::Index numVecs = std::min<Eigen::Index>(dim, requestedVecs);
    if (numVecs <= 0) return 0.0;
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(dim, numVecs);
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(Q);
    Q = qr.householderQ() * Eigen::MatrixXd::Identity(dim, numVecs);

    Eigen::VectorXd prevEvals = Eigen::VectorXd::Zero(numVecs);
    
    for (int iter = 0; iter < 200; ++iter) { 
        Eigen::MatrixXd V(dim, numVecs);

        #pragma omp parallel for schedule(static)
        for (int col = 0; col < numVecs; ++col) {
            V.col(col) = solver.solve(Q.col(col));
        }
        Eigen::HouseholderQR<Eigen::MatrixXd> qrIter(V);
        Q = qrIter.householderQ() * Eigen::MatrixXd::Identity(dim, numVecs);
        
        Eigen::MatrixXd SQ_iter = S * Q;
        Eigen::VectorXd currentEvals = (Q.transpose() * SQ_iter).diagonal();
        
        if (iter > 0 && (currentEvals - prevEvals).norm() < 1e-6) {
            break; 
        }
        prevEvals = currentEvals;
    }

    Eigen::MatrixXd SQ = S * Q;
    Eigen::MatrixXd T = Q.transpose() * SQ;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(T);
    if (es.info() != Eigen::Success) {
        return 0.0;
    }

    Eigen::VectorXd evals = es.eigenvalues();
    if (evals.size() == 0) return 0.0;

    Eigen::Index idx = std::min(modeIndex, static_cast<std::size_t>(evals.size() - 1));
    double lambda = evals(idx) ;
    return std::sqrt(std::max(0.0, lambda));
}
