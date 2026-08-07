#include "massMatrix.h"

Eigen::SparseMatrix<double> getSparseMassMatrix(
    std::size_t dim,
    const double uniteMass,
    const double radius,
    const int massNum
) {
    const int totalMassNum{massNum * massNum * massNum};
    Eigen::SparseMatrix<double> M(dim, dim);
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(dim);

    for (int i{0}; i < totalMassNum; ++i) {
        const int base{6 * i};
        triplets.emplace_back(base + 0, base + 0, uniteMass);
        triplets.emplace_back(base + 1, base + 1, uniteMass);
        triplets.emplace_back(base + 2, base + 2, uniteMass);
        triplets.emplace_back(base + 3, base + 3, uniteMass);
        const double rotInertia = uniteMass * (2.0 / 5.0) * radius * radius;
        triplets.emplace_back(base + 4, base + 4, rotInertia);
        triplets.emplace_back(base + 5, base + 5, rotInertia);
    }
    M.setFromTriplets(triplets.begin(), triplets.end());
    return M;
}