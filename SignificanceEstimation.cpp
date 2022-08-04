#include "SignificanceEstimation.h"


Alignment::Alignment() {
    alignmentFileName = "";
    threshold = 0;
    numberOfStrings = 0;
    lengthOfAlignment = 0;
    lengthOfSeedAlignment = 0;
};

Alignment::Alignment(std::string alignmentFileName, double threshold) :
        alignmentFileName(alignmentFileName), threshold(threshold) {
    std::ifstream file("../" + alignmentFileName);
    assert(file.is_open());

    std::string tmpString;
    file >> tmpString;
    assert(!file.eof());
    file.seekg(0);
    lengthOfAlignment = tmpString.length();
    columns.resize(lengthOfAlignment);
    std::vector<int> countOfGaps(lengthOfAlignment);
    numberOfStrings = 0;
    do {
        file >> tmpString;
        ++numberOfStrings;
        for(int i = 0; i < lengthOfAlignment; ++i)
            countOfGaps[i] = tmpString[i] == '-' ? countOfGaps[i] + 1 : countOfGaps[i];

    } while(!file.eof());
    lengthOfSeedAlignment = 0;
    for(int i = 0; i < lengthOfAlignment; ++i) {
        columns[i] = (double) countOfGaps[i] / numberOfStrings < threshold;    // false = ignored column.
        lengthOfSeedAlignment = columns[i] ? ++lengthOfSeedAlignment : lengthOfSeedAlignment;
    }

    file.close();
}
std::string Alignment::getAlignmentFileName() {
    return alignmentFileName;
}

double Alignment::getThreshold() {
    return threshold;
}

int Alignment::getNumberOfStrings() {
    return numberOfStrings;
}

int Alignment::getLengthOfAlignment() {
    return lengthOfAlignment;
}

int Alignment::getLengthOfSeedAlignment() {
    return lengthOfSeedAlignment;
}

std::vector<bool> Alignment::getColumns() {
    return columns;
}
//----------------------------------------------------------------------------------------------------------------------
std::map<char, double> BackgroundModel::getEmissions() {
    return emissions;
}
//----------------------------------------------------------------------------------------------------------------------
void PHMM::PHMMCount(Alignment alignment) {
    std::ifstream file("../" + alignment.getAlignmentFileName());
    assert(file.is_open());

    int lengthOfSeedAlignment = alignment.getLengthOfSeedAlignment();
    std::string tmpString;
    std::vector<bool> columns = alignment.getColumns();
    int lengthOfAlignment = alignment.getLengthOfAlignment();
    do {
        file >> tmpString;
        char tmpTypeOfState = 'S';
        for(int stringIndex = 0, PHMMIndex = 0; stringIndex < lengthOfAlignment; ++stringIndex) {
            // Transitions counting.
            if(!columns[stringIndex]) {
                if(tmpString[stringIndex] != '-') {
                    switch(tmpTypeOfState) {
                        case 'I':
                            ++transitions[PHMMIndex][PHMMIndex];
                            break;
                        case 'M':
                            ++transitions[PHMMIndex][PHMMIndex + 2];
                            PHMMIndex += 2;
                            break;
                        default:
                            ++transitions[PHMMIndex][PHMMIndex + 1];
                            PHMMIndex += 1;
                    }
                    tmpTypeOfState = 'I';
                }
            }
            else {
                if(tmpString[stringIndex] == '-') {
                    switch(tmpTypeOfState) {
                        case 'I':
                            ++transitions[PHMMIndex][PHMMIndex + 2];
                            PHMMIndex += 2;
                            break;
                        case 'M':
                            ++transitions[PHMMIndex][PHMMIndex + 4];
                            PHMMIndex += 4;
                            break;
                        default:
                            ++transitions[PHMMIndex][PHMMIndex + 3];
                            PHMMIndex += 3;
                    }
                    tmpTypeOfState = 'D';
                }
                else {
                    switch(tmpTypeOfState) {
                        case 'I':
                            ++transitions[PHMMIndex][PHMMIndex + 1];
                            PHMMIndex += 1;
                            break;
                        case 'M':
                            ++transitions[PHMMIndex][PHMMIndex + 3];
                            PHMMIndex += 3;
                            break;
                        default:
                            ++transitions[PHMMIndex][PHMMIndex + 2];
                            PHMMIndex += 2;
                    }
                    tmpTypeOfState = 'M';
                }
            }
            // Emissions counting.
            if(tmpString[stringIndex] != '-')
                ++emissions[PHMMIndex][tmpString[stringIndex]];
        }
        // The path ends in the terminal state.
        if(tmpTypeOfState == 'M')
            ++transitions[lengthOfSeedAlignment * 3 - 1][lengthOfSeedAlignment * 3 + 2];
        else if(tmpTypeOfState == 'D')
            ++transitions[lengthOfSeedAlignment * 3][lengthOfSeedAlignment * 3 + 2];
        else if(tmpTypeOfState == 'I')
            ++transitions[lengthOfSeedAlignment * 3 + 1][lengthOfSeedAlignment * 3 + 2];
    } while(!file.eof());

    file.close();
}

void PHMM::PHMMPseudocountsAndNormalizationTr(Alignment alignment, double pseudocountValue) {
    double sum;
    // Add pseudocounts to the first rectangle.
    int xStart = 1, yStart = 0;
    for(int y = yStart; y < 2; ++y) {
        sum = 0;
        for (int x = xStart; x < 4; ++x) {
            if (transitions[y][x] == 0)
                transitions[y][x] = pseudocountValue;
            sum += transitions[y][x];
        }
        for (int x = xStart; x < 4; ++x)
            transitions[y][x] /= sum;
    }
    // Add pseudocounts to the middle rectangles.
    xStart += 3; yStart += 2;
    int lengthOfSeedAlignment = alignment.getLengthOfSeedAlignment();
    while(xStart <= lengthOfSeedAlignment * 3 - 1) {
        for(int y = yStart; y < yStart + 3; ++y) {
            sum = 0;
            for (int x = xStart; x < xStart + 3; ++x) {
                if (transitions[y][x] == 0)
                    transitions[y][x] = pseudocountValue;
                sum += transitions[y][x];
            }
            for (int x = xStart; x < xStart + 3; ++x)
                transitions[y][x] /= sum;
        }
        xStart += 3; yStart += 3;
    }
    // Add pseudocounts to the last rectangle.
    for(int y = yStart; y < yStart + 3; ++y) {
        sum = 0;
        for (int x = xStart; x < xStart + 2; ++x) {
            if (transitions[y][x] == 0)
                transitions[y][x] = pseudocountValue;
            sum += transitions[y][x];
        }
        for (int x = xStart; x < xStart + 2; ++x)
            transitions[y][x] /= sum;
    }
}

void PHMM::PHMMPseudocountsAndNormalizationEm(Alignment alignment, double pseudocountValue) {
    double sum;
    int lengthOfSeedAlignment = alignment.getLengthOfSeedAlignment();
    for(int y = 1; y < lengthOfSeedAlignment * 3 + 2; ++y) {
        if(y % 3 == 0) {
            emissions[y]['-'] = 1;
            continue;
        }
        sum = 0;
        for (auto& x: emissions[y]) {
            if (x.first != '-' && x.second == 0)
                x.second = pseudocountValue;
            sum += x.second;
        }
        for (auto& x: emissions[y])
            x.second /= sum;
    }
}

PHMM::PHMM(Alignment alignment, double pseudocountValue) {
    // Creating tables.
    int lengthOfSeedAlignment = alignment.getLengthOfSeedAlignment();
    transitions.resize(lengthOfSeedAlignment * 3 + 3);
    emissions.resize(lengthOfSeedAlignment * 3 + 3);
    for(int i = 0; i < lengthOfSeedAlignment * 3 + 3; ++i) {
        transitions[i].resize(lengthOfSeedAlignment * 3 + 3);
        emissions[i]['-'] = 0; emissions[i]['A'] = 0; emissions[i]['C'] = 0;
        emissions[i]['D'] = 0; emissions[i]['E'] = 0; emissions[i]['F'] = 0;
    }

    // Count transitions and emissions.
    PHMMCount(alignment);
    // Add pseudocounts and normalize for transitions.
    PHMMPseudocountsAndNormalizationTr(alignment, pseudocountValue);
    // Add pseudocounts and normalize for emissions.
    PHMMPseudocountsAndNormalizationEm(alignment, pseudocountValue);
}

std::vector<std::vector<double>> PHMM::getTransitions() {
    return transitions;
}

std::vector<std::map<char, double>> PHMM::getEmissions() {
    return emissions;
}
//----------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------
SignificanceEstimation::SignificanceEstimation(std::string alignmentFileName, double threshold, double pseudocountValue) :
        alignment(alignmentFileName, threshold), phmm(alignment, pseudocountValue) {

}

Alignment SignificanceEstimation::getAlignment() {
    return alignment;
}

BackgroundModel SignificanceEstimation::getBackgroundModel() {
    return backgroundModel;
}

PHMM SignificanceEstimation::getPhmm() {
    return phmm;
}

Sample SignificanceEstimation::getSample() {
    return sample;
}

std::vector<std::vector<double>> SignificanceEstimation::getTransitionsForSample() {
    return transitionsForSample;
}

std::vector<std::map<char, double>> SignificanceEstimation::getEmissionsForSample() {
    return emissionsForSample;
}

double  SignificanceEstimation::partitionFunction(std::string sequence, double T) {
    assert(T > 0); double Tln = 1 / T;
    int lengthOfSeedAlignment = alignment.getLengthOfSeedAlignment(),
            lengthOfSequence = alignment.getLengthOfAlignment();
    std::vector<std::vector<double>> transitions = phmm.getTransitions();
    std::vector<std::map<char, double>> emissions = phmm.getEmissions();
    std::vector<std::vector<double>> forward(lengthOfSeedAlignment * 3 + 3);
    for(std::vector<double>& v: forward)
        v.resize(lengthOfSequence);

    // Initialisation. Filling the first 3 lines.
    forward[1][0] = pow(transitions[0][1], Tln) * pow(emissions[1][sequence[0]], Tln);
    forward[2][0] = pow(transitions[0][2], Tln) * pow(emissions[2][sequence[0]], Tln);
    forward[3][0] = pow(transitions[0][3], Tln) * pow(emissions[3][sequence[0]], Tln);
    for(int x = 1; x < lengthOfSequence; ++x){
        forward[1][x] = forward[1][x - 1] * pow(transitions[1][1], Tln) * pow(emissions[1][sequence[x]], Tln);
        forward[2][x] = forward[1][x - 1] * pow(transitions[1][2], Tln) * pow(emissions[2][sequence[x]], Tln);
        forward[3][x] = forward[1][x - 1] * pow(transitions[1][3], Tln) * pow(emissions[3][sequence[x]], Tln);
    }

    // Middle calculations.
    for(int y = 4; y < lengthOfSeedAlignment * 3 + 1; y += 3) {
        for(int x = 1; x < lengthOfSequence; ++x) {
            forward[y][x] = (forward[y - 2][x - 1] * pow(transitions[y - 2][y], Tln) +
                    forward[y - 1][x - 1] * pow(transitions[y - 1][y], Tln) +
                    forward[y][x - 1] * pow(transitions[y][y], Tln)) *
                    pow(emissions[y][sequence[x]], Tln);
            forward[y + 1][x] = (forward[y - 2][x - 1] * pow(transitions[y - 2][y + 1], Tln) +
                    forward[y - 1][x - 1] * pow(transitions[y - 1][y + 1], Tln) +
                    forward[y][x - 1] * pow(transitions[y][y + 1], Tln)) *
                    pow(emissions[y + 1][sequence[x]], Tln);
            forward[y + 2][x] = (forward[y - 2][x - 1] * pow(transitions[y - 2][y + 2], Tln) +
                    forward[y - 1][x - 1] * pow(transitions[y - 1][y + 2], Tln) +
                    forward[y][x - 1] * pow(transitions[y][y + 2], Tln)) *
                    pow(emissions[y + 2][sequence[x]], Tln);
        }
    }

    // Completion.
    for(int x = 1; x < lengthOfSequence; ++x) {
        forward[lengthOfSeedAlignment * 3 + 1][x] = (forward[lengthOfSeedAlignment * 3 - 1][x - 1] *
                pow(transitions[lengthOfSeedAlignment * 3 - 1][lengthOfSeedAlignment - 2], Tln) +
                forward[lengthOfSeedAlignment * 3][x - 1] *
                pow(transitions[lengthOfSeedAlignment * 3][lengthOfSeedAlignment - 2], Tln) +
                forward[lengthOfSeedAlignment * 3 + 1][x - 1] *
                pow(transitions[lengthOfSeedAlignment * 3 + 1][lengthOfSeedAlignment - 2], Tln)) *
                pow(emissions[lengthOfSeedAlignment * 3 + 1][sequence[x]], Tln);
    }

//    for(auto f: forward) {
//        for(auto inf: f){
//            std::cout << inf << ' ';
//        }
//        std::cout << std::endl;
//    }

    return forward[lengthOfSeedAlignment * 3 - 1][lengthOfSequence - 1] +
           forward[lengthOfSeedAlignment * 3][lengthOfSequence - 1] +
           forward[lengthOfSeedAlignment * 3 + 1][lengthOfSequence - 1];
}

void SignificanceEstimation::averageEmissionsCalculation(double T) {
    assert(T > 0); double Tln = 1 / T;
    std::map<char, double> emissionsBM = backgroundModel.getEmissions();
    int lengthOfSeedAlignment = alignment.getLengthOfSeedAlignment();
    std::vector<std::map<char, double>> emissions = phmm.getEmissions();
    averageEmissions.resize(lengthOfSeedAlignment * 3 + 3);

    for(int y = 0; y < lengthOfSeedAlignment * 3 + 3; ++y) {
        averageEmissions[y] = emissionsBM['-'] * pow(emissions[y]['-'], Tln) +
                emissionsBM['A'] * pow(emissions[y]['A'], Tln) +
                emissionsBM['C'] * pow(emissions[y]['C'], Tln) +
                emissionsBM['D'] * pow(emissions[y]['D'], Tln) +
                emissionsBM['E'] * pow(emissions[y]['E'], Tln) +
                emissionsBM['F'] * pow(emissions[y]['F'], Tln);
    }
}

double SignificanceEstimation::ZCalculation(double T) {
    assert(T > 0); double Tln = 1 / T;
    int lengthOfSeedAlignment = alignment.getLengthOfSeedAlignment(),
            lengthOfSequence = alignment.getLengthOfAlignment();
    averageEmissionsCalculation(T);
    std::vector<std::vector<double>> transitions = phmm.getTransitions();
    transitionsForSample.resize(lengthOfSeedAlignment * 3 + 3);
    for(std::vector<double>& v: transitionsForSample)
        v.resize(lengthOfSequence);

    // Initialisation. Filling the first 3 lines.
    transitionsForSample[1][0] = pow(transitions[0][1], Tln) * averageEmissions[1];
    transitionsForSample[2][0] = pow(transitions[0][2], Tln) * averageEmissions[2];
    transitionsForSample[3][0] = pow(transitions[0][3], Tln) * averageEmissions[3];
    for(int x = 1; x < lengthOfSequence; ++x){
        transitionsForSample[1][x] = transitionsForSample[1][x - 1] * pow(transitions[1][1], Tln) * averageEmissions[1];
        transitionsForSample[2][x] = transitionsForSample[1][x - 1] * pow(transitions[1][2], Tln) * averageEmissions[2];
        transitionsForSample[3][x] = transitionsForSample[1][x - 1] * pow(transitions[1][3], Tln) * averageEmissions[3];
    }

    // Middle calculations.
    for(int y = 4; y < lengthOfSeedAlignment * 3 + 1; y += 3) {
        for(int x = 1; x < lengthOfSequence; ++x) {
            transitionsForSample[y][x] = (transitionsForSample[y - 2][x - 1] * pow(transitions[y - 2][y], Tln) +
                    transitionsForSample[y - 1][x - 1] * pow(transitions[y - 1][y], Tln) +
                    transitionsForSample[y][x - 1] * pow(transitions[y][y], Tln)) *
                    averageEmissions[y];
            transitionsForSample[y + 1][x] = (transitionsForSample[y - 2][x - 1] * pow(transitions[y - 2][y + 1], Tln) +
                    transitionsForSample[y - 1][x - 1] * pow(transitions[y - 1][y + 1], Tln) +
                    transitionsForSample[y][x - 1] * pow(transitions[y][y + 1], Tln)) *
                    averageEmissions[y + 1];
            transitionsForSample[y + 2][x] = (transitionsForSample[y - 2][x - 1] * pow(transitions[y - 2][y + 2], Tln) +
                    transitionsForSample[y - 1][x - 1] * pow(transitions[y - 1][y + 2], Tln) +
                    transitionsForSample[y][x - 1] * pow(transitions[y][y + 2], Tln)) *
                    averageEmissions[y + 2];
        }
    }

    // Completion.
    for(int x = 1; x < lengthOfSequence; ++x) {
        transitionsForSample[lengthOfSeedAlignment * 3 + 1][x] = (transitionsForSample[lengthOfSeedAlignment * 3 - 1][x - 1] *
                pow(transitions[lengthOfSeedAlignment * 3 - 1][lengthOfSeedAlignment - 2], Tln) +
                transitionsForSample[lengthOfSeedAlignment * 3][x - 1] *
                pow(transitions[lengthOfSeedAlignment * 3][lengthOfSeedAlignment - 2], Tln) +
                transitionsForSample[lengthOfSeedAlignment * 3 + 1][x - 1] *
                pow(transitions[lengthOfSeedAlignment * 3 + 1][lengthOfSeedAlignment - 2], Tln)) *
                averageEmissions[lengthOfSeedAlignment * 3 + 1];
    }

//    for(auto f: transitionsForSample) {
//        for(auto inf: f){
//            std::cout << inf << ' ';
//        }
//        std::cout << std::endl;
//    }

    Z = transitionsForSample[lengthOfSeedAlignment * 3 - 1][lengthOfSequence - 1] +
        transitionsForSample[lengthOfSeedAlignment * 3][lengthOfSequence - 1] +
        transitionsForSample[lengthOfSeedAlignment * 3 + 1][lengthOfSequence - 1];
    return Z;
}

void SignificanceEstimation::emissionsForSampleCalculation(double T) {
    assert(T > 0); double Tln = 1 / T;
    int lengthOfSeedAlignment = alignment.getLengthOfSeedAlignment();
    std::map<char, double> emissionsBM = backgroundModel.getEmissions();
    std::vector<std::map<char, double>> emissions = phmm.getEmissions();
    emissionsForSample.resize(lengthOfSeedAlignment);

    for(int i = 0; i < lengthOfSeedAlignment * 3 + 3; ++i) {
        emissionsForSample[i]['-'] = emissionsBM['-'] * pow(emissions[i]['-'], Tln) / averageEmissions[i];
        emissionsForSample[i]['A'] = emissionsBM['A'] * pow(emissions[i]['A'], Tln) / averageEmissions[i];
        emissionsForSample[i]['C'] = emissionsBM['C'] * pow(emissions[i]['C'], Tln) / averageEmissions[i];
        emissionsForSample[i]['D'] = emissionsBM['D'] * pow(emissions[i]['D'], Tln) / averageEmissions[i];
        emissionsForSample[i]['E'] = emissionsBM['E'] * pow(emissions[i]['E'], Tln) / averageEmissions[i];
        emissionsForSample[i]['F'] = emissionsBM['F'] * pow(emissions[i]['F'], Tln) / averageEmissions[i];
    }
}