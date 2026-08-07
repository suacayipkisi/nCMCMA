cmake_minimum_required(VERSION 3.22)
project(ncmcma VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

set(SOURCES
	main.cpp
	src/matrixOperations/stdEigenValueSolver.cpp
	src/modalAnalysis/massMatrix.cpp
	src/modalAnalysis/stiffMatrix.cpp
)

set(HEADERS
	src/modalAnalysis/massMatrix.h
	src/modalAnalysis/stiffMatrix.h
	src/matrixOperations/stdEigenValueSolver.h
)

add_compile_definitions(MAIN_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

add_compile_options(-O3 -march=native)
add_link_options(-O3)

add_executable(ncmcma ${SOURCES})
target_link_libraries(ncmcma PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX)
