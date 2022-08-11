#include "SignificanceEstimation.h"


int main() {
    std::cout.setf(std::ios::fixed);
    double T = 1;
    std::cout.precision(10);
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
    std::cout << std::endl;
    std::cout << std::endl;
    for(auto i: emissions) {
        for (auto j: i) {
            std::cout << std::setw(4) <<  j.second << ' ';
        }
        std::cout << std::endl;
    }
//    BackgroundModel backgroundModel = significanceEstimation.getBackgroundModel();
//    std::cout << std::endl << significanceEstimation.partitionFunction("ADDEFAAADF", T) << std::endl;
//    std::cout << std::endl << backgroundModel.probabilityOfString("ADDEFAAADF") << std::endl;
    std::cout << std::endl << significanceEstimation.ZCalculation(10, T) << std::endl;
    significanceEstimation.emissionsForSampleCalculation(T);
    Sample sample = significanceEstimation.getSample();
    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();
    sample.sampleEmission(2, emissionsForSample);
    sample.sampleSequences(10, 10, 8,
                           emissionsForSample, transitionsForSample);


//    std::cout << std::endl << significanceEstimation.ZCalculation(T) << std::endl;
//    significanceEstimation.emissionsForSampleCalculation(T);
//    Sample sample = significanceEstimation.getSample();
//    Alignment alignment = significanceEstimation.getAlignment();
//    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();
//    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
//    for(int i = 0; i < 80; ++i) {
//        sample.sampleSequence(alignment, emissionsForSample, transitionsForSample);
//    }


    return 0;
};
