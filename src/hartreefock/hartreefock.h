#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include <armadillo>
#include <mpi.h>
#include "math.h"
#include "system/system.h"

using namespace arma;

#define RUN_MPI

class HartreeFock
{
public:
    HartreeFock(System *system);

    virtual void solve()=0;
    double getEnergy();
    field<mat> getQmatrix();
    virtual field<mat> getCoeff()=0;
    virtual field<colvec> getFockEnergy()=0;

protected:
    System *m_system;

    mat m_h;
    mat m_S;
    field<mat> m_Q;
    mat m_V;

    int m_matDim;
    double m_energy;

    double m_toler;
    double m_restrictedFactor;

    virtual void buildFockMatrix()=0;
    void calcIntegrals();
    void diagOverlap();
    void solveSingle(const mat &Fock, mat &Coeffs, mat &P, colvec &fockEnergy, int nElectrons);
};

#endif // HARTREEFOCK_H
