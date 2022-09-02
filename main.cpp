#include "SignificanceEstimation.h"


int main() {
//    std::cout.setf(std::ios::fixed);
//    std::cout.precision(10);
    double T = 4, threshold = pow(10, -12);
    int lengthOfSequence = 8, numberOfSequences = 1000;
    SignificanceEstimation significanceEstimation("alignmentFile.txt", "sampleFile.txt",
                                                  0.4, 0.01);

    double Z = significanceEstimation.ZCalculation(lengthOfSequence, T);
    significanceEstimation.emissionsForSampleCalculation(T);
    Sample &sample = significanceEstimation.getSample();
    Alignment alignment = significanceEstimation.getAlignment();
    std::vector<std::vector<double>> transitionsForSample = significanceEstimation.getTransitionsForSample();
    std::vector<std::map<char, double>> emissionsForSample = significanceEstimation.getEmissionsForSample();

    sample.allSequencesGeneration(8);
    std::cout << "fpr: " << significanceEstimation.fprCalculation(threshold) << std::endl;


//    std::cout << "T: " << significanceEstimation.temperatureChoice(lengthOfSequence, threshold) << std::endl;
    int count = 10;
    double sumfpr = 0, tmp1;
    for(int i = 0; i < count; ++i) {
        sample.sampleSequences(numberOfSequences, lengthOfSequence, alignment.getLengthOfSeedAlignment(),
                               emissionsForSample, transitionsForSample);
        std::cout << "fprEstimation: " << (tmp1 = significanceEstimation.fprEstimation(threshold, Z, T)) << std::endl;
        sumfpr += tmp1;
    }
    std::cout << std::endl << "average fprEstimation: " << sumfpr / count << std::endl;







    return 0;
};
