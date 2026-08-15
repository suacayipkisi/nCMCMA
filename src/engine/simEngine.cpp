#include "simEngine.h"

#include <iostream>
#include <cmath>
#include <numbers>
#include <mutex>

#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#else
#error "Eigen headers not found. Add Eigen include path (e.g., /usr/include/eigen3)."
#endif

#include "../modalAnalysis/stiffMatrix.h"
#include "../modalAnalysis/massMatrix.h"
#include "../matrixOperations/stdEigenValueSolver.h"
#include "../parameters.h"

void runSimulation(SimulationParameters& g_params, SimulationState& g_simState) {
    g_simState.isRunning = true;
    g_simState.progress = 0.0f;

    {
        std::lock_guard<std::mutex> lock(g_simState.dataMutex);
        g_simState.results.clear();
    }

    const double mass = (4.0 / 3.0) * std::numbers::pi * std::pow(g_params.radius, 3) * g_params.density;
    const std::size_t dim = 6 * g_params.massNum * g_params.massNum * g_params.massNum;

    std::cout << "[ENGINE] Creating sparse mass and stiffness matrix...\n";
    Eigen::SparseMatrix<double> sparseM = getSparseMassMatrix(dim, mass, g_params.radius, g_params.massNum);
    Eigen::SparseMatrix<double> sparseK = getSparseStiffnessMatrix(dim, g_params.stiffnessConst, g_params.radius, g_params.massNum);

    int shift = 0;
    bool equalSituation = true;
    double wOld = -0.002;
    int iter = 0;

    std::cout << "[ENGINE] Frequency searching starting. Aim: " << g_params.maxTargetHz << " Hz\n";

    while (true) {
        while (equalSituation) {
            const double wNow = getNaturalFrequency(sparseM, sparseK, (iter + shift));
            if (std::abs(wOld - wNow) < 0.001) {
                shift++;
            } else {
                equalSituation = false;
                wOld = wNow;
            }
        }
        equalSituation = true;

        const double w = wOld;
        const double currentHz = w / (2.0 * std::numbers::pi);

        if (currentHz > g_params.maxTargetHz) {
            std::cout << "[ENGINE] Reached maximum frequency: " << currentHz << " Hz. Calculation stopping.\n";
            break;
        }

        std::cout << "[ENGINE] Finded frequency: " << currentHz << " Hz (" << w << " rad/s)\n";

        const double w2 = w * w;
        std::complex<double> k_coeff(1.0, w * g_params.beta);
        std::complex<double> m_coeff(-w2, w * g_params.alpha);

        Eigen::SparseMatrix<std::complex<double>> Z = 
            sparseK.cast<std::complex<double>>() * k_coeff + 
            sparseM.cast<std::complex<double>>() * m_coeff;

        Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>> lu;
        lu.compute(Z);
        
        if (lu.info() != Eigen::Success) {
            std::cerr << "[ERROR] Z matrix SparseLU factorization is unsuccessful!\n";
            break;
        }

        Eigen::VectorXcd Q_vec = Eigen::VectorXcd::Zero(dim);
        Q_vec(0) = std::complex<double>(1.0, 0.0);
        Eigen::VectorXcd q_vec = lu.solve(Q_vec);

        {
            std::lock_guard<std::mutex> lock(g_simState.dataMutex);
            g_simState.results.push_back({
                currentHz, 
                std::vector<std::complex<double>>(q_vec.begin(), q_vec.end())
            });
        }

        iter++;
        
        g_simState.progress = static_cast<float>(currentHz / g_params.maxTargetHz);
    }

    g_simState.isRunning = false;
    std::cout << "[ENGINE] Simulation done.\n";
}
