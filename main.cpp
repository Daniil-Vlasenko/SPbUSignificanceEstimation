#include "SignificanceEstimation.h"


int main() {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(10);
    SignificanceEstimation significanceEstimation("alignmentFile.txt", 0.4, 0.05);
    std::vector<std::vector<double>> transitions = significanceEstimation.getPhmm().getTransitions();
    std::vector<std::map<char, double>> emissions = significanceEstimation.getPhmm().getEmissions();
    for(auto i: transitions) {
        for (auto j: i) {
            std::cout << std::setw(4) << j << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    for(auto i: emissions) {
        for (auto j: i) {
            std::cout << std::setw(4) <<  j.second << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << significanceEstimation.partitionFunction("ACDEFACADF", 1) << std::endl;
    std::cout << std::endl;
    std::cout << significanceEstimation.ZCalculation(1) << std::endl;

    return 0;
};
