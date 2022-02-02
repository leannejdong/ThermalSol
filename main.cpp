#include <iostream>
#include "ThermalSolve.h"
int main() {
    try {
        ThermalSol::ThermalSolve();
    } catch (const std::exception &exception) {
        std::cerr << exception.what() << "\n";
        throw;
    }
    return 0;
}
