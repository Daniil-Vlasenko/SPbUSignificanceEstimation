#include "SignificanceEstimation.h"


int main() {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(10);
    double T = 1, threshold = pow(10, -3), significanceLevel = 0.99;
    int lengthOfSequence = 8, numberOfSequences = 10000;
    SignificanceEstimation significanceEstimation("alignmentFile.txt", "sampleFile.txt",
                                                  0.4, 0.01);


    double Z = significanceEstimation.ZCalculation(lengthOfSequence, T);
    significanceEstimation.emissionsForSampleCalculation(T);
    Sample &sample = significanceEstimation.getSample();
    Alignment alignment = significanceEstimation.getAlignment();
    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();
    sample.sampleSequences(numberOfSequences, lengthOfSequence, alignment.getLengthOfSeedAlignment(),
                           emissionsForSample, transitionsForSample);
    std::pair<double, double> confidenceInterval =
            significanceEstimation.confidenceIntervalCalculation(threshold, Z, T, significanceLevel);
    std::cout << confidenceInterval.first << " " << confidenceInterval.second << std::endl;
    double fpr1 = significanceEstimation.fprEstimation(threshold, Z, T);
//    sample.allSequencesGeneration(8);
//    double fpr2 = significanceEstimation.fprCalculation(threshold);


//    meanTheta in fprCalculation: 0.00013312







    return 0;
};
