#include "stiffMatrix.h"

void fillDiag(
    const std::array<int, 3>& elementPosition,
    const int elementID,
    std::vector<std::vector<double>>& stiffnessMatrix,
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
    switch (situation)
    {
    case 0:
        stiffnessMatrix[base + 0][base + 0] = 2 * stiffConst;
        stiffnessMatrix[base + 1][base + 1] = 2 * stiffConst;
        stiffnessMatrix[base + 2][base + 2] = 2 * stiffConst;
        stiffnessMatrix[base + 3][base + 3] = 4 * stiffConst * radius * radius;
        stiffnessMatrix[base + 4][base + 4] = 4 * stiffConst * radius * radius;
        stiffnessMatrix[base + 5][base + 5] = 4 * stiffConst * radius * radius;
        break;
    case 1:
        stiffnessMatrix[base + 0][base + 0] = 2 * stiffConst;
        stiffnessMatrix[base + 1][base + 1] = 2 * stiffConst;
        stiffnessMatrix[base + 2][base + 2] = 2 * stiffConst;
        stiffnessMatrix[base + 3][base + 3] = 4 * stiffConst * radius * radius;
        stiffnessMatrix[base + 4][base + 4] = 3 * stiffConst * radius * radius;
        stiffnessMatrix[base + 5][base + 5] = 3 * stiffConst * radius * radius;
        break;
    case 2:
        stiffnessMatrix[base + 0][base + 0] = 1 * stiffConst;
        stiffnessMatrix[base + 1][base + 1] = 1 * stiffConst;
        stiffnessMatrix[base + 2][base + 2] = 2 * stiffConst;
        stiffnessMatrix[base + 3][base + 3] = 3 * stiffConst * radius * radius;
        stiffnessMatrix[base + 4][base + 4] = 3 * stiffConst * radius * radius;
        stiffnessMatrix[base + 5][base + 5] = 2 * stiffConst * radius * radius;
        break;
    case 3:
        stiffnessMatrix[base + 0][base + 0] = stiffConst;
        stiffnessMatrix[base + 1][base + 1] = stiffConst;
        stiffnessMatrix[base + 2][base + 2] = stiffConst;
        stiffnessMatrix[base + 3][base + 3] = 2 * stiffConst * radius * radius;
        stiffnessMatrix[base + 4][base + 4] = 2 * stiffConst * radius * radius;
        stiffnessMatrix[base + 5][base + 5] = 2 * stiffConst * radius * radius;
        break;
    default:
        break;
    }
}

void fillElement(
    const std::array<int, 3>& elementPosition,
    const int elementID,
    const std::array<int, 3>& effectorPosition,
    const int effectorID,
    std::vector<std::vector<double>>& stiffnessMatrix,
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
        stiffnessMatrix[nPi + 0][nPj + 0] = (stiffConst * (-1));
        stiffnessMatrix[nPi + 1][nPj + 1] = (stiffConst * (-1));
        stiffnessMatrix[nPi + 2][nPj + 2] = (stiffConst * (-1));
        stiffnessMatrix[nPi + 4][nPj + 4] = stiffConst * (radius * radius * (-1));
        stiffnessMatrix[nPi + 5][nPj + 5] = stiffConst * (radius * radius *(-1));
        if(state == 1){
            stiffnessMatrix[nPi + 1][nPj + 5] = stiffConst * (radius);
            stiffnessMatrix[nPi + 2][nPj + 4] = stiffConst * (radius * (-1));
            stiffnessMatrix[nPi + 4][nPj + 2] = stiffConst * (radius * (-1));
            stiffnessMatrix[nPi + 5][nPj + 1] = stiffConst * (radius);
        }
        else if(state == -1){
            stiffnessMatrix[nPi + 1][nPj + 5] = stiffConst * (radius * (-1));
            stiffnessMatrix[nPi + 2][nPj + 4] = stiffConst * (radius);
            stiffnessMatrix[nPi + 4][nPj + 2] = stiffConst * (radius);
            stiffnessMatrix[nPi + 5][nPj + 1] = stiffConst * (radius * (-1));
        }
    }
    else if(
        elementPosition[0] == effectorPosition[0]
        && elementPosition[2] == effectorPosition[2]
    ){
        int state{effectorPosition[1] - elementPosition[1]};
        stiffnessMatrix[nPi + 0][nPj + 0] = (stiffConst * (-1));
        stiffnessMatrix[nPi + 1][nPj + 1] = (stiffConst * (-1));
        stiffnessMatrix[nPi + 2][nPj + 2] = (stiffConst * (-1));
        stiffnessMatrix[nPi + 3][nPj + 3] = stiffConst * (radius * radius * (-1));
        stiffnessMatrix[nPi + 5][nPj + 5] = stiffConst * (radius * radius *(-1));
        if(state == 1){
            stiffnessMatrix[nPi + 0][nPj + 5] = stiffConst * (radius * (-1));
            stiffnessMatrix[nPi + 2][nPj + 3] = stiffConst * (radius);
            stiffnessMatrix[nPi + 5][nPj + 0] = stiffConst * (radius * (-1));
            stiffnessMatrix[nPi + 3][nPj + 2] = stiffConst * (radius);
        }
        else if(state == -1){
            stiffnessMatrix[nPi + 2][nPj + 3] = stiffConst * (radius * (-1));
            stiffnessMatrix[nPi + 5][nPj + 0] = stiffConst * (radius);
            stiffnessMatrix[nPi + 0][nPj + 5] = stiffConst * (radius);
            stiffnessMatrix[nPi + 3][nPj + 2] = stiffConst * (radius * (-1));
        }
    }
    else{
        int state{effectorPosition[2] - elementPosition[2]};
        stiffnessMatrix[nPi + 0][nPj + 0] = (stiffConst * (-1));
        stiffnessMatrix[nPi + 1][nPj + 1] = (stiffConst * (-1));
        stiffnessMatrix[nPi + 2][nPj + 2] = (stiffConst * (-1));
        stiffnessMatrix[nPi + 3][nPj + 3] = stiffConst * (radius * radius * (-1));
        stiffnessMatrix[nPi + 4][nPj + 4] = stiffConst * (radius * radius *(-1));
        if(state == 1){
            stiffnessMatrix[nPi + 1][nPj + 3] = stiffConst * (radius * (-1));
            stiffnessMatrix[nPi + 0][nPj + 4] = stiffConst * (radius);
            stiffnessMatrix[nPi + 3][nPj + 1] = stiffConst * (radius * (-1));
            stiffnessMatrix[nPi + 4][nPj + 0] = stiffConst * (radius);
        }
        else if(state == -1){
            stiffnessMatrix[nPi + 0][nPj + 4] = stiffConst * (radius * (-1));
            stiffnessMatrix[nPi + 1][nPj + 3] = stiffConst * (radius);
            stiffnessMatrix[nPi + 3][nPj + 1] = stiffConst * (radius);
            stiffnessMatrix[nPi + 4][nPj + 0] = stiffConst * (radius * (-1));
        }
    }

    return;
}

std::vector<std::vector<double>> getStiffnessMatrix(
    std::size_t dim,
    const double stiffConst,
    const double radius,
    const int massNum)
{
    const int totalMassNum{massNum * massNum * massNum};
    std::vector<std::vector<double>> stiffnessMatrix(dim, std::vector<double>(dim, 0.0));
    std::array<int, 3> elePos{0, 0, 0};
    std::array<int, 3> effPos{0, 0, 0};

    const std::array<std::array<int, 3>, 6> offsetDirs = {{
        {1, 0, 0}, {-1, 0, 0},
        {0, 1, 0}, {0, -1, 0},
        {0, 0, 1}, {0, 0, -1}
    }};

    for (int eleID{1}; eleID <= totalMassNum; ++eleID){
        int idx = eleID - 1;
        int N2 = massNum * massNum;
        elePos[0] = idx % massNum;
        elePos[1] = (idx % N2) / massNum;
        elePos[2] = idx / N2;
        fillDiag(elePos, eleID, stiffnessMatrix, stiffConst, radius, massNum);

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
                fillElement(elePos, eleID, effPos, effID, stiffnessMatrix, stiffConst, radius);
            }
        }
    }

    return stiffnessMatrix;
}

Eigen::SparseMatrix<double> getSparseStiffnessMatrix(
    std::size_t dim,
    const double stiffConst,
    const double radius,
    const int massNum
) {
    const int totalMassNum{massNum * massNum * massNum};
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(dim * 45);

    std::vector<std::vector<double>> denseK = getStiffnessMatrix(dim, stiffConst, radius, massNum);

    for (std::size_t i = 0; i < dim; ++i) {
        for (std::size_t j = 0; j < dim; ++j) {
            if (denseK[i][j] != 0.0) {
                triplets.emplace_back(i, j, denseK[i][j]);
            }
        }
    }

    Eigen::SparseMatrix<double> K(dim, dim);
    K.setFromTriplets(triplets.begin(), triplets.end());
    return K;
}