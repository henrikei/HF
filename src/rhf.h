#ifndef RHF_H
#define RHF_H

#include "hartreefock.h"

class RHF : public HartreeFock
{
public:
    RHF(System *newSystem);
    void solve();
    mat getCoeff();
private:
    mat F;
    mat C;
    mat P;
    colvec fockEnergy;

    void buildMatrix();
};

#endif // RHF_H
