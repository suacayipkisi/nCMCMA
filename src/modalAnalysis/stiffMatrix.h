#pragma once

#include <array>
#include <cstddef>

// real mass num is pow(massNum, 3)


template <std::size_t T>
void fillDiag(
    const std::array<int, 3>& elementPosition,
    const int elementID,
    std::array<std::array<double, T>, T>& stiffnessMatrix,
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

template <std::size_t T>
void fillElement(
    const std::array<int, 3>& elementPosition,
    const int elementID,
    const std::array<int, 3>& effectorPosition,
    const int effectorID,
    std::array<std::array<double, T>, T>& stiffnessMatrix,
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

template <std::size_t T>
std::array<std::array<double, T>, T> getStriffnessMatrix(
    const double stiffConst,
    const double radius,
    const int massNum)
{
    const int totalMassNum{massNum * massNum * massNum};
    static std::array<std::array<double, T>, T> stiffnessMatrix{};
    std::array<int, 3> elePos{0, 0, 0};
    std::array<int, 3> effPos{0, 0, 0};
    for (int eleID{1}; eleID <= totalMassNum; ++eleID){
        int idx = eleID - 1;
        int N2 = massNum * massNum;
        elePos[0] = idx % massNum;
        elePos[1] = (idx % N2) / massNum;
        elePos[2] = idx / N2;
        fillDiag(elePos, eleID, stiffnessMatrix, stiffConst, radius, massNum);
        for (int effID{1}; effID <= totalMassNum; ++effID){
            int fidx = effID - 1;
            int fN2 = massNum * massNum;
            effPos[0] = fidx % massNum;
            effPos[1] = (fidx % fN2) / massNum;
            effPos[2] = fidx / fN2;
            if(eleID != effID){
                fillElement(elePos, eleID, effPos, effID, stiffnessMatrix, stiffConst, radius);
            }
        }
    }

    static std::array<std::array<double, T>, T> arrStiffnessMatrix{};
    for (std::size_t i{0}; i < static_cast<std::size_t>(T); ++i){
        for (std::size_t j{0}; j < static_cast<std::size_t>(T); ++j){
            arrStiffnessMatrix[i][j] = stiffnessMatrix[i][j];
        }
    }
    return arrStiffnessMatrix;
}