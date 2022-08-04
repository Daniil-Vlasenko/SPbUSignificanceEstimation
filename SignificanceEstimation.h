#pragma once
#include <iostream>
#include<iomanip>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <cmath>


class Alignment {
private:
    std::string alignmentFileName;
    double threshold;
    int numberOfStrings;
    int lengthOfAlignment;
    int lengthOfSeedAlignment;
    std::vector<bool> columns;    // false = ignored column.
public:
    // Identification of ignored strings
    Alignment();
    Alignment(std::string alignmentFileName, double threshold);
    std::string getAlignmentFileName();
    int getNumberOfStrings();
    int getLengthOfAlignment();
    int getLengthOfSeedAlignment();
    std::vector<bool> getColumns();
};

class BackgroundModel {
private:
    std::map<char, double> emissions = {{'-', 1./6}, {'A', 1./6}, {'C', 1./6},
                                        {'D', 1./6}, {'E', 1./6}, {'F', 1./6},};
public:
    std::map<char, double> getEmissions();
};

class PHMM {
private:
    std::vector<std::vector<double>> transitions;
    std::vector<std::map<char, double>> emissions;

    // Called in the constructor to count transitions and emissions.
    void PHMMCount(Alignment alignment);
    // Called in the constructor to add pseudocounts and normalize counts for transitions.
    void PHMMPseudocountsAndNormalizationTr(Alignment alignment, double pseudocountValue);
    // Called in the constructor to add pseudocounts and normalize counts for emissions.
    void PHMMPseudocountsAndNormalizationEm(Alignment alignment, double pseudocountValue);

public:
    PHMM(Alignment alignment, double pseudocountValue);
    std::vector<std::vector<double>> getTransitions();
    std::vector<std::map<char, double>> getEmissions();
};

class Sample {
private:
    std::vector<std::vector<double>> transitions;
public:
    std::vector<std::vector<double>> getTransitions();
};

class SignificanceEstimation {
    Alignment alignment;
    BackgroundModel backgroundModel;
    PHMM phmm;
    Sample sample;

};
