#include "SignificanceEstimation.h"


int main() {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(10);
    SignificanceEstimation significanceEstimation("alignmentFile.txt", "sampleFile.txt", 0.4, 0.05);
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
    std::cout << std::endl;
    Sample sample = significanceEstimation.getSample();
    significanceEstimation.emissionsForSampleCalculation(1);
    Alignment alignment = significanceEstimation.getAlignment();
    std::cout << sample.sampleSequence(alignment,
                                       significanceEstimation.getEmissionsForSample(),
                                       significanceEstimation.getTransitionsForSample()) << std::endl;
    sample.sampleSequences(10, alignment,
           significanceEstimation.getEmissionsForSample(),
           significanceEstimation.getTransitionsForSample());

    return 0;
};
