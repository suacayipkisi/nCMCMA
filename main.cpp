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
    const int numThreads{static_cast<int>(hardwareThreads)};
    omp_set_num_threads(numThreads/2);
    Eigen::setNbThreads(numThreads/2);

    SimulationParameters parameters{};

    SimulationState state{};
    
    std::cout << "Thread number in system: " << (hardwareThreads)
              << " | Using thread number: " << numThreads << "\n\n";

    initAndRunGui(parameters, state);

    return 0;
}
