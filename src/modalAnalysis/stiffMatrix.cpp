#include <omp.h>

#include "stiffMatrix.h"

static void fillDiagTriplets(
    const std::array<int, 3>& elementPosition,
    const int elementID,
    std::vector<Eigen::Triplet<double>>& triplets,
    const double stiffConst,
    const double radius,
    const int massNum
){
    const int base{6 * (elementID - 1)};

    int situation{};
    for (int i{0}; i < 3; ++i){
        if((elementPosition[i] == 0) || (elementPosition[i] == massNum - 1)){
            situation +=1;
        }
    }
    const double r2 = radius * radius;
    switch (situation)
    {
    case 0:
        triplets.emplace_back(base + 0, base + 0, 2 * stiffConst);
        triplets.emplace_back(base + 1, base + 1, 2 * stiffConst);
        triplets.emplace_back(base + 2, base + 2, 2 * stiffConst);
        triplets.emplace_back(base + 3, base + 3, 4 * stiffConst * r2);
        triplets.emplace_back(base + 4, base + 4, 4 * stiffConst * r2);
        triplets.emplace_back(base + 5, base + 5, 4 * stiffConst * r2);
        break;
    case 1:
        triplets.emplace_back(base + 0, base + 0, 2 * stiffConst);
        triplets.emplace_back(base + 1, base + 1, 2 * stiffConst);
        triplets.emplace_back(base + 2, base + 2, 2 * stiffConst);
        triplets.emplace_back(base + 3, base + 3, 4 * stiffConst * r2);
        triplets.emplace_back(base + 4, base + 4, 3 * stiffConst * r2);
        triplets.emplace_back(base + 5, base + 5, 3 * stiffConst * r2);
        break;
    case 2:
        triplets.emplace_back(base + 0, base + 0, 1 * stiffConst);
        triplets.emplace_back(base + 1, base + 1, 1 * stiffConst);
        triplets.emplace_back(base + 2, base + 2, 2 * stiffConst);
        triplets.emplace_back(base + 3, base + 3, 3 * stiffConst * r2);
        triplets.emplace_back(base + 4, base + 4, 3 * stiffConst * r2);
        triplets.emplace_back(base + 5, base + 5, 2 * stiffConst * r2);
        break;
    case 3:
        triplets.emplace_back(base + 0, base + 0, stiffConst);
        triplets.emplace_back(base + 1, base + 1, stiffConst);
        triplets.emplace_back(base + 2, base + 2, stiffConst);
        triplets.emplace_back(base + 3, base + 3, 2 * stiffConst * r2);
        triplets.emplace_back(base + 4, base + 4, 2 * stiffConst * r2);
        triplets.emplace_back(base + 5, base + 5, 2 * stiffConst * r2);
        break;
    default:
        break;
    }
}

static void fillElementTriplets(
    const std::array<int, 3>& elementPosition,
    const int elementID,
    const std::array<int, 3>& effectorPosition,
    const int effectorID,
    std::vector<Eigen::Triplet<double>>& triplets,
    const double stiffConst,
    const double radius
)
{
    const int nPi{6 * (elementID - 1)};
    const int nPj{6 * (effectorID - 1)};

    if(
        elementPosition[1] == effectorPosition[1]
        && elementPosition[2] == effectorPosition[2]
    ){
        int state{effectorPosition[0] - elementPosition[0]};
        triplets.emplace_back(nPi + 0, nPj + 0, -stiffConst);
        triplets.emplace_back(nPi + 1, nPj + 1, -stiffConst);
        triplets.emplace_back(nPi + 2, nPj + 2, -stiffConst);
        triplets.emplace_back(nPi + 4, nPj + 4, -stiffConst * radius * radius);
        triplets.emplace_back(nPi + 5, nPj + 5, -stiffConst * radius * radius);
        if(state == 1){
            triplets.emplace_back(nPi + 1, nPj + 5, stiffConst * radius);
            triplets.emplace_back(nPi + 2, nPj + 4, -stiffConst * radius);
            triplets.emplace_back(nPi + 4, nPj + 2, -stiffConst * radius);
            triplets.emplace_back(nPi + 5, nPj + 1, stiffConst * radius);
        }
        else if(state == -1){
            triplets.emplace_back(nPi + 1, nPj + 5, -stiffConst * radius);
            triplets.emplace_back(nPi + 2, nPj + 4, stiffConst * radius);
            triplets.emplace_back(nPi + 4, nPj + 2, stiffConst * radius);
            triplets.emplace_back(nPi + 5, nPj + 1, -stiffConst * radius);
        }
    }
    else if(
        elementPosition[0] == effectorPosition[0]
        && elementPosition[2] == effectorPosition[2]
    ){
        int state{effectorPosition[1] - elementPosition[1]};
        triplets.emplace_back(nPi + 0, nPj + 0, -stiffConst);
        triplets.emplace_back(nPi + 1, nPj + 1, -stiffConst);
        triplets.emplace_back(nPi + 2, nPj + 2, -stiffConst);
        triplets.emplace_back(nPi + 3, nPj + 3, -stiffConst * radius * radius);
        triplets.emplace_back(nPi + 5, nPj + 5, -stiffConst * radius * radius);
        if(state == 1){
            triplets.emplace_back(nPi + 0, nPj + 5, -stiffConst * radius);
            triplets.emplace_back(nPi + 2, nPj + 3, stiffConst * radius);
            triplets.emplace_back(nPi + 5, nPj + 0, -stiffConst * radius);
            triplets.emplace_back(nPi + 3, nPj + 2, stiffConst * radius);
        }
        else if(state == -1){
            triplets.emplace_back(nPi + 2, nPj + 3, -stiffConst * radius);
            triplets.emplace_back(nPi + 5, nPj + 0, stiffConst * radius);
            triplets.emplace_back(nPi + 0, nPj + 5, stiffConst * radius);
            triplets.emplace_back(nPi + 3, nPj + 2, -stiffConst * radius);
        }
    }
    else{
        int state{effectorPosition[2] - elementPosition[2]};
        triplets.emplace_back(nPi + 0, nPj + 0, -stiffConst);
        triplets.emplace_back(nPi + 1, nPj + 1, -stiffConst);
        triplets.emplace_back(nPi + 2, nPj + 2, -stiffConst);
        triplets.emplace_back(nPi + 3, nPj + 3, -stiffConst * radius * radius);
        triplets.emplace_back(nPi + 4, nPj + 4, -stiffConst * radius * radius);
        if(state == 1){
            triplets.emplace_back(nPi + 1, nPj + 3, -stiffConst * radius);
            triplets.emplace_back(nPi + 0, nPj + 4, stiffConst * radius);
            triplets.emplace_back(nPi + 3, nPj + 1, -stiffConst * radius);
            triplets.emplace_back(nPi + 4, nPj + 0, stiffConst * radius);
        }
        else if(state == -1){
            triplets.emplace_back(nPi + 0, nPj + 4, -stiffConst * radius);
            triplets.emplace_back(nPi + 1, nPj + 3, stiffConst * radius);
            triplets.emplace_back(nPi + 3, nPj + 1, stiffConst * radius);
            triplets.emplace_back(nPi + 4, nPj + 0, -stiffConst * radius);
        }
    }
}

Eigen::SparseMatrix<double> getSparseStiffnessMatrix(
    std::size_t dim,
    const double stiffConst,
    const double radius,
    const int massNum,
    const int numThreads
)
{
    const int totalMassNum{massNum * massNum * massNum};
    const int N2 = massNum * massNum;

    std::vector<std::vector<Eigen::Triplet<double>>> threadTriplets(numThreads);

    const std::array<std::array<int, 3>, 6> offsetDirs = {{
        {1, 0, 0}, {-1, 0, 0},
        {0, 1, 0}, {0, -1, 0},
        {0, 0, 1}, {0, 0, -1}
    }};

    #pragma omp parallel num_threads(numThreads)
    {
        const int threadID = omp_get_thread_num();
        auto& localTriplets = threadTriplets[threadID];
        localTriplets.reserve((dim * 45) / numThreads);

        std::array<int, 3> elePos{0, 0, 0};
        std::array<int, 3> effPos{0, 0, 0};

        #pragma omp for schedule(static)
        for (int eleID = 1; eleID <= totalMassNum; ++eleID) {
            int idx = eleID - 1;
            elePos[0] = idx % massNum;
            elePos[1] = (idx % N2) / massNum;
            elePos[2] = idx / N2;
            fillDiagTriplets(elePos, eleID, localTriplets, stiffConst, radius, massNum);

            for (const auto& dir : offsetDirs) {
                effPos[0] = elePos[0] + dir[0];
                effPos[1] = elePos[1] + dir[1];
                effPos[2] = elePos[2] + dir[2];

                if (effPos[0] >= 0 && effPos[0] < massNum &&
                    effPos[1] >= 0 && effPos[1] < massNum &&
                    effPos[2] >= 0 && effPos[2] < massNum)
                {
                    int fidx = effPos[0] + effPos[1] * massNum + effPos[2] * N2;
                    int effID = fidx + 1;
                    fillElementTriplets(elePos, eleID, effPos, effID, localTriplets, stiffConst, radius);
                }
            }
        }
    }
    
    std::size_t totalTripletsSize = 0;
    for (const auto& tVec : threadTriplets) {
        totalTripletsSize += tVec.size();
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(totalTripletsSize);

    for (auto& tVec : threadTriplets) {
        triplets.insert(triplets.end(), 
                        std::make_move_iterator(tVec.begin()), 
                        std::make_move_iterator(tVec.end()));
    }

    Eigen::SparseMatrix<double> K(dim, dim);
    K.setFromTriplets(triplets.begin(), triplets.end());
    K.makeCompressed();
    return K;
}
