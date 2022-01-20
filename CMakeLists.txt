cmake_minimum_required(VERSION 3.22)
project(ThermalSol)

set(CMAKE_CXX_STANDARD 20)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(ThermalSol main.cpp)
target_link_libraries(ThermalSol Eigen3::Eigen)
