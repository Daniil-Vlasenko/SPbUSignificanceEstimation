#include "SignificanceEstimation.h"


int main() {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(10);
    double T = 10;
    int lengthOfSequence = 8, numberOfSequences = 1000;
    SignificanceEstimation significanceEstimation("alignmentFile.txt", "sampleFile.txt",
                                                  0.4, 0.01);
//    std::vector<std::vector<double>> transitions = significanceEstimation.getPhmm().getTransitions();
//    std::vector<std::map<char, double>> emissions = significanceEstimation.getPhmm().getEmissions();
//    for(auto i: transitions) {
//        for (auto j: i) {
//            std::cout << std::setw(4) << j << ' ';
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl << std::endl;
//    for(auto i: emissions) {
//        for (auto j: i) {
//            std::cout << std::setw(4) <<  j.second << ' ';
//        }
//        std::cout << std::endl;
//    }



    double Z = significanceEstimation.ZCalculation(lengthOfSequence, T);
    significanceEstimation.emissionsForSampleCalculation(T);
    Sample &sample = significanceEstimation.getSample();
    Alignment alignment = significanceEstimation.getAlignment();
    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();

//    double fprSum = 0;
    int count = 10;
    for(int i = 0; i < count; ++i) {
        sample.sampleSequences(numberOfSequences, lengthOfSequence, alignment.getLengthOfSeedAlignment(),
                               emissionsForSample, transitionsForSample);
        std::cout << significanceEstimation.fprCalculation(pow(10, -5), Z, T) << std::endl;
    }
//    std::cout << fprSum / count;





    return 0;
};
