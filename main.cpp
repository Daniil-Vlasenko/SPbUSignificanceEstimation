#include "SignificanceEstimation.h"


int main() {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(10);
    double T = 1, threshold = pow(10, -85), significanceLevel = 0.99;
    int lengthOfSequence = 100, numberOfSequences = 10000;
    SignificanceEstimation significanceEstimation("alignmentFile.txt", "sampleFile.txt",
                                                  0.4, 0.01);
    Sample &sample = significanceEstimation.getSample();
    Alignment alignment = significanceEstimation.getAlignment();
    BackgroundModel backgroundModel = significanceEstimation.getBackgroundModel();

    double Z = significanceEstimation.ZCalculation(lengthOfSequence, T);
    significanceEstimation.emissionsForSampleCalculation(T);
    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();
    //---
    sample.sampleSequencesB(numberOfSequences, lengthOfSequence, backgroundModel);
    std::pair<double, double> confidenceInterval =
            significanceEstimation.confidenceIntervalCalculation(threshold, significanceLevel);
    std::cout << confidenceInterval.first << " " << confidenceInterval.second << std::endl;
    //---
//    std::cout << "T: " << significanceEstimation.temperatureChoice(lengthOfSequence, threshold) << std::endl;
    sample.sampleSequences(numberOfSequences, lengthOfSequence, alignment.getLengthOfSeedAlignment(),
                           emissionsForSample, transitionsForSample);
    double fpr1 = significanceEstimation.fprEstimation(threshold, Z, T);
    std::cout << fpr1 << std::endl;
    //---
//    sample.allSequencesGeneration(8);
//    double fpr2 = significanceEstimation.fprCalculation(threshold);
    //---










//for(int t = 98; t < 105; ++t) {
//    threshold = pow(10, -t);
//    sample.sampleSequencesB(numberOfSequences, lengthOfSequence, backgroundModel);
//    std::pair<double, double> confidenceInterval =
//            significanceEstimation.confidenceIntervalCalculation(threshold, significanceLevel);
//    std::cout << "10e-" << t << ": [" << confidenceInterval.first << "; " << confidenceInterval.second << "]" << std::endl;
//}








    return 0;
};
