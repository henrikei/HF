#ifndef RHF_H
#define RHF_H

#include "hartreefock.h"

class RHF : public HartreeFock
{
public:
    RHF(System *system);
    void solve();
    field<mat> getCoeff();
    field<mat> getDensityMatrix();
    field<colvec> getFockEnergy();

private:
    mat m_F;
    mat m_C;
    mat m_P;
    colvec m_fockEnergy;
    int m_nElectrons;

    void buildFockMatrix();
};

#endif // RHF_H
