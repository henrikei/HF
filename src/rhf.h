#ifndef RHF_H
#define RHF_H

#include "hartreefock.h"

class RHF : public HartreeFock
{
public:
    RHF(System *system, int perturbOrder=1);
    void solve();
    mat getCoeff();

private:
    mat m_F;
    mat m_C;
    mat m_P;
    colvec m_fockEnergy;
    int m_nElectrons;

    field<mat> m_MOI;

    void buildFockMatrix();
    double perturbation2order();
};

#endif // RHF_H
