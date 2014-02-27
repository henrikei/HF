#ifndef UHF_H
#define UHF_H

#include "hartreefock.h"


class UHF : public HartreeFock
{
public:
    UHF(System *system);
    UHF(System *system, int perturbOrder);
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

    void buildFockMatrix();
    double perturbation2order();
};

#endif // UHF_H
