#pragma once
#include <iostream>
#include<iomanip>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <cmath>
#include <chrono>
#include <random>
#include <algorithm>


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
    double getThreshold();
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
    double probabilityOfString(std::string sequence);
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
    unsigned seed;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
    std::string sampleFileName;
    int numberOfSequences;

public:
    Sample(std::string sampleFileName);
    std::string getSampleFileName();
    int getNumberOfSequences();
    // Sample Emission for state(I, M, D).
    char sampleEmission(int state, const std::vector<std::map<char, double>> &emissionsForSample);
    std::string sampleSequence(Alignment alignment, const std::vector<std::map<char, double>> &emissionsForSample,
                               const std::vector<std::vector<double>> &transitionsForSample);
    std::string sampleSequences(int numberOfSequences, Alignment alignment,
                                const std::vector<std::map<char, double>> &emissionsForSample,
                                const std::vector<std::vector<double>> &transitionsForSample);

};

class SignificanceEstimation {
private:
    Alignment alignment;
    BackgroundModel backgroundModel;
    PHMM phmm;
    Sample sample;
    std::vector<double> averageEmissions;
    std::vector<std::vector<double>> transitionsForSample;
    std::vector<std::map<char, double>> emissionsForSample;

public:
    SignificanceEstimation(std::string alignmentFileName, std::string sampleFileName, double threshold, double pseudocountValue);
    Alignment getAlignment();
    BackgroundModel getBackgroundModel();
    PHMM getPhmm();
    Sample& getSample();
    std::vector<std::vector<double>> getTransitionsForSample();
    std::vector<std::map<char, double>> getEmissionsForSample();
    // Calculation of Z(D,T). If T = 1, result is the probabilityOfString of sequence.
    double partitionFunction(std::string sequence, double T);
    // Calculation of average emissions for Z(T).
    void averageEmissionsCalculation(double T);
    // Calculation of Z(T) and transitionsForSample.
    double ZCalculation(double T);
    // Calculation of emissionsForSample.
    void emissionsForSampleCalculation(double T);
    // Calculation of false positive rate.
    double fprCalculation(double threshold, double Z);

};
