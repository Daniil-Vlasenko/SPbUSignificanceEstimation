#include "SignificanceEstimation.h"


int main() {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(10);
    double T = 2;
    int lengthOfSequence = 10, numberOfSequences = 100;
    SignificanceEstimation significanceEstimation("alignmentFile.txt", "sampleFile.txt",
                                                  0.4, 0.005);
    std::vector<std::vector<double>> transitions = significanceEstimation.getPhmm().getTransitions();
    std::vector<std::map<char, double>> emissions = significanceEstimation.getPhmm().getEmissions();
    for(auto i: transitions) {
        for (auto j: i) {
            std::cout << std::setw(4) << j << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
    for(auto i: emissions) {
        for (auto j: i) {
            std::cout << std::setw(4) <<  j.second << ' ';
        }
        std::cout << std::endl;
    }

    std::cout << std::endl << significanceEstimation.partitionFunction("AFDEFADADC", T) << std::endl << std::endl;
    std::cout << significanceEstimation.ZCalculation(10, T) << std::endl << std::endl;
    significanceEstimation.emissionsForSampleCalculation(1);
    Sample &sample = significanceEstimation.getSample();
    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();
    sample.sampleSequences(100, 10, 8, emissionsForSample, transitionsForSample);

//    double Z = significanceEstimation.ZCalculation(lengthOfSequence, T);
//    significanceEstimation.emissionsForSampleCalculation(T);
//    Sample &sample = significanceEstimation.getSample();
//    Alignment alignment = significanceEstimation.getAlignment();
//    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
//    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();
//
//    for(int i = 0; i < 10; ++i) {
//        sample.sampleSequences(numberOfSequences, lengthOfSequence, alignment.getLengthOfSeedAlignment(),
//                               emissionsForSample, transitionsForSample);
//        std::cout << significanceEstimation.fprCalculation(0.00001, Z, T) << std::endl;
//    }





    return 0;
};
