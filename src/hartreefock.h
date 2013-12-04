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
    virtual void solve()=0;
    double getEnergy();
    virtual mat getCoeff()=0;
protected:
    System *system;

    mat h;
    double**** Q;
    mat S;
    int matDim;
    int nElectrons;
    double energy;
    double toler;

    virtual void buildMatrix()=0;
    void calcIntegrals();
    void solveSingle(const mat &Fock, mat &Coeffs, mat &P, colvec &fockEnergy);
};

#endif // HARTREEFOCK_H
