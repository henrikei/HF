#ifndef UHF_H
#define UHF_H

#include "hartreefock.h"


class UHF : public HartreeFock
{
public:
    UHF(System *system, int perturbOrder=1);
    void solve();
    mat getCoeff();

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

    field<mat> m_MOI_UU;
    field<mat> m_MOI_DD;
    field<mat> m_MOI_UD;

    void buildFockMatrix();
    double perturbation2order();
};

#endif // UHF_H
