#ifndef UHF_H
#define UHF_H

#include "hartreefock.h"


class UHF : public HartreeFock
{
public:
    UHF(System *newSystem);
    void solve();
    mat getCoeff();
private:
    mat Fup;
    mat Fdown;
    mat Cup;
    mat Cdown;
    mat Pup;
    mat Pdown;
    colvec fockEnergyUp;
    colvec fockEnergyDown;

    void buildMatrix();
};

#endif // UHF_H
