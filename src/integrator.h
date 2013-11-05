#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <armadillo>
#include <vector>
#include "boysfunction.h"

using namespace arma;
using namespace std;



class Integrator
{
public:
    Integrator(int angMomMax);
    ~Integrator();

    void setAlpha(double newAlpha);
    void setBeta(double newBeta);
    void setGamma(double newGamma);
    void setDelta(double newDelta);
    void setPositionA(rowvec3 RA);
    void setPositionB(rowvec3 RB);
    void setPositionC(rowvec3 RC);
    void setPositionD(rowvec3 RD);
    void setMaxAngMom(int newAngMom);
    void setE_AB();
    void setE_CD();
    void setR(double a, rowvec3 A, int coulombType);
    double overlap(int i, int j, int k, int l, int m, int n);
    double kinetic(int i, int j, int k, int l, int m, int n);
    double coulomb1(int i, int j, int k, int l, int m, int n);
    double coulomb2(int i1, int j1, int k1, int l1, int m1, int n1, int i2, int j2, int k2, int l2, int m2, int n2);

    cube E_AB[3];
    cube E_CD[3];
    vector<cube> R;

private:
    BoysFunction *boys;
    mat Rnuclei;
    double alpha, beta, gamma, delta;
    double angMom;
};

#endif // INTEGRATOR_H
