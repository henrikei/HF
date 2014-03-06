#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include <armadillo>
#include "math.h"
#include "system.h"

using namespace arma;


class HartreeFock
{
public:
    HartreeFock(System *system);
    virtual void solve()=0;
    double getEnergy();
    double getEnergyMP2();
    virtual mat getCoeff()=0;
protected:
    System *m_system;
    mat m_h;
    mat m_S;
    field<mat> m_Q;
    int m_matDim;
    double m_energy;
    double m_energyMP2;
    double m_toler;
    double m_restrictedFactor;
    int m_perturbOrder;

    virtual void buildFockMatrix()=0;
    void calcIntegrals();
    void solveSingle(const mat &Fock, mat &Coeffs, mat &P, colvec &fockEnergy, int nElectrons);
    virtual double perturbation2order()=0;
};

#endif // HARTREEFOCK_H
