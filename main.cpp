#include "SignificanceEstimation.h"


int main() {
    std::cout.setf(std::ios::fixed);
    double T = 1;
    std::cout.precision(10);
    SignificanceEstimation significanceEstimation("alignmentFile.txt", "sampleFile.txt",
                                                  0.4, 0.05);
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
    std::cout << std::endl << significanceEstimation.partitionFunction("AC-FEFADCF", T) << std::endl;
    std::cout << std::endl << significanceEstimation.ZCalculation(T) << std::endl;
    significanceEstimation.emissionsForSampleCalculation(T);
    Sample sample = significanceEstimation.getSample();
    Alignment alignment = significanceEstimation.getAlignment();
    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();
    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
    for(int i = 0; i < 80; ++i) {
        sample.sampleSequence(alignment, emissionsForSample, transitionsForSample);
    }




//    std::cout << std::endl;
//    std::cout << significanceEstimation.partitionFunction("ACDEFACADF", 1) << std::endl;
//    std::cout << std::endl;
//    std::cout << significanceEstimation.ZCalculation(1) << std::endl;
//    std::cout << std::endl;
//    significanceEstimation.ZCalculation(1);
//    significanceEstimation.emissionsForSampleCalculation(1);
//    Alignment alignment = significanceEstimation.getAlignment();
//    significanceEstimation.getSample().sampleSequences(10, alignment,
//           significanceEstimation.getEmissionsForSample(),
//           significanceEstimation.getTransitionsForSample());
//    std::cout << significanceEstimation.fprCalculation(0.0000000000000005);

    return 0;
};
