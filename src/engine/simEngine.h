#pragma once

#include <atomic>
#include <mutex>
#include <vector>
#include <complex>

struct SimulationParameters {
    int massNum{4};
    float alpha{0.02287f};
    float beta{0.00000578f};
    float radius{0.01f};
    float density{7800.0f};
    float stiffnessConst{1000.0f};
    float maxTargetHz{50.0f};

    SimulationParameters() = default;
};

struct ModeResult {
    double hz;
    std::vector<std::complex<double>> displacements;
};

struct SimulationState {
    std::atomic<bool> isRunning{false};
    std::atomic<float> progress{0.0f};
    std::mutex dataMutex;
    
    std::vector<ModeResult> results; 

    SimulationState() = default;
    ~SimulationState() = default;
    SimulationState(const SimulationState&) = delete;
    SimulationState& operator=(const SimulationState&) = delete;
    SimulationState(SimulationState&&) = delete;
    SimulationState& operator=(SimulationState&&) = delete;
};

extern SimulationParameters g_params;
extern SimulationState g_simState;

void runSimulation();