#pragma once

#include <atomic>
#include <complex>
#include <mutex>
#include <vector>

struct ModeResult {
    double hz;
    std::vector<std::complex<double>> displacements;
};

struct SimulationParameters {
    int massNum{};
    float alpha{};
    float beta{};
    float radius{};
    float density{};
    float stiffnessConst{};
    float maxTargetHz{};

    SimulationParameters():
        massNum{4},
        alpha{0.02287f},
        beta{0.00000578f},
        radius{0.01f},
        density{7800.0f},
        stiffnessConst{1000.0f},
        maxTargetHz{50.0f}
    {}
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
