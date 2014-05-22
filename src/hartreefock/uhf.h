#ifndef UHF_H
#define UHF_H

#include "hartreefock.h"


class UHF : public HartreeFock
{
public:
    UHF(System *system);
    void solve();
    field<mat> getCoeff();
    field<mat> getDensityMatrix();
    field<colvec> getFockEnergy();
    int getNumOfElectronsUp();
    int getNumOfElecrtonsDown();
    double getSpinExpectation();

private:
    mat m_Fup;
    mat m_Fdown;
    mat m_Cup;
    mat m_Cdown;
    mat m_Pup;
    mat m_Pdown;
    colvec m_fockEnergyUp;
    colvec m_fockEnergyDown;
    int m_nElectronsUp;
    int m_nElectronsDown;

    void buildFockMatrix();
};

#endif // UHF_H
