#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include <armadillo>
#include "math.h"
#include "system.h"

using namespace arma;


class HartreeFock
{
public:
    HartreeFock(System *newSystem);
    void solve();
    double getEnergy();
    vec getCoeff();
private:
    System *system;

    mat h;
    double**** Q;
    mat F;
    mat S;
    vec C;
    int matDim;
    double fockEnergy;
    double energy;
    double toler;

    void buildMatrix();
    void calcIntegrals();
    void solveSingle();
};

#endif // HARTREEFOCK_H
