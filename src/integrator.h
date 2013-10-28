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
    Integrator();
    ~Integrator();

    void setAlpha(double newAlpha);
    void setBeta(double newBeta);
    void setPositionA(rowvec3 RA);
    void setPositionB(rowvec3 RB);
    void setMaxAngMom(int newAngMom);
    void setE();
    void setR(rowvec3 C);
    double overlap(int i, int j, int k, int l, int m, int n);
    double kinetic(int i, int j, int k, int l, int m, int n);
    double coulomb1(int i, int j, int k, int l, int m, int n);

    cube E[3];
    vector<cube> R;

private:
    mat Rnuclei;
    double alpha, beta;
    double angMom;
};

#endif // INTEGRATOR_H
