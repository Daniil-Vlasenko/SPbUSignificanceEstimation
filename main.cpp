#include "SignificanceEstimation.h"


int main() {
//    std::cout.setf(std::ios::fixed);
//    std::cout.precision(10);
    double T = 1;
    int lengthOfSequence = 10, numberOfSequences = 10000;
    SignificanceEstimation significanceEstimation("alignmentFile.txt", "sampleFile.txt",
                                                  0.4, 0.005);
//    std::vector<std::vector<double>> transitions = significanceEstimation.getPhmm().getTransitions();
//    std::vector<std::map<char, double>> emissions = significanceEstimation.getPhmm().getEmissions();
//    for(auto i: transitions) {
//        for (auto j: i) {
//            std::cout << std::setw(4) << j << ' ';
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;
//    std::cout << std::endl;
//    for(auto i: emissions) {
//        for (auto j: i) {
//            std::cout << std::setw(4) <<  j.second << ' ';
//        }
//        std::cout << std::endl;
//    }
//    BackgroundModel backgroundModel = significanceEstimation.getBackgroundModel();
//    std::cout << std::endl << significanceEstimation.partitionFunction("ADDEFAAADF", T) << std::endl;
//    std::cout << std::endl << backgroundModel.probabilityOfString("ADDEFAAADF") << std::endl;
    double Z = significanceEstimation.ZCalculation(lengthOfSequence, T);
    significanceEstimation.emissionsForSampleCalculation(T);
    Sample &sample = significanceEstimation.getSample();
    Alignment alignment = significanceEstimation.getAlignment();
    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();


    for(int i = 0; i < 10; ++i) {
        sample.sampleSequences(numberOfSequences, lengthOfSequence, alignment.getLengthOfSeedAlignment(),
                               emissionsForSample, transitionsForSample);
        std::cout << significanceEstimation.fprCalculation(pow(10, -100), Z, T) << std::endl;
    }





    return 0;
};
