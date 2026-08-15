#include <omp.h>
#include <thread>
#include <iostream>

#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#elif __has_include(<eigen3/Eigen/Dense>)
#include <eigen3/Eigen/Dense>
#endif

#include "src/parameters.h"
#include "src/gui/gui.h"

int main() {
    const unsigned int hardwareThreads{std::thread::hardware_concurrency()};
    const int activeThreads = std::max(1, static_cast<int>(hardwareThreads) / 2);
    
    omp_set_num_threads(activeThreads);
    Eigen::setNbThreads(activeThreads);

    SimulationParameters parameters{
        activeThreads, 
        4, 
        0.02287f, 
        0.00000578f, 
        0.01f, 
        800.0f, 
        1000.0f, 
        50.0f
    };

    SimulationState state{};
    
    std::cout << "Thread number in system: " << (hardwareThreads)
              << " | Using thread number: " << activeThreads << "\n\n";

    initAndRunGui(parameters, state);

    return 0;
}
