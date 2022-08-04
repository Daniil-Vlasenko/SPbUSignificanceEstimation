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
std::vector<std::vector<double>> Sample::getTransitions() {
    return transitions;
}
//----------------------------------------------------------------------------------------------------------------------