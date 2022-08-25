#include "SignificanceEstimation.h"


int main() {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(10);
    double T = 1, threshold = pow(10, -3);
    int lengthOfSequence = 10, numberOfSequences = 1000;
    SignificanceEstimation significanceEstimation("alignmentFile.txt", "sampleFile.txt",
                                                  0.4, 0.01);

    double Z = significanceEstimation.ZCalculation(lengthOfSequence, T);
    significanceEstimation.emissionsForSampleCalculation(T);
    Sample &sample = significanceEstimation.getSample();
    Alignment alignment = significanceEstimation.getAlignment();
    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();


    std::cout << significanceEstimation.temperatureChoice(lengthOfSequence, threshold) << std::endl;
    int count = 10;
    double sum = 0, tmp;
    for(int i = 0; i < count; ++i) {
        sample.sampleSequences(numberOfSequences, lengthOfSequence, alignment.getLengthOfSeedAlignment(),
                               emissionsForSample, transitionsForSample);
        std::cout << (tmp = significanceEstimation.fprCalculation(threshold, Z, T)) << std::endl;
        sum += tmp;
    }
    std::cout << std::endl << sum / count;







    return 0;
};
