#pragma once

class Probdist {
public:
    double pNormal(double x);
    double xNormal(double prob);

    double pChi(double x, int n);
    double xChi(double prob, int n);
private:
    const double Eps = 1e-15;

    int fequal( double a, double b );
    int fcompare( double a, double b );
    void NORMAL( int type, double &x, double &p);
    void  CHI( int type, double n, double &x, double &p);

};
