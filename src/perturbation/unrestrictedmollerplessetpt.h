#ifndef UNRESTRICTEDMOLLERPLESSETPT_H
#define UNRESTRICTEDMOLLERPLESSETPT_H

#include <armadillo>
#include "mollerplessetpt.h"
#include "uhf.h"

using namespace std;
using namespace arma;


class UnrestrictedMollerPlessetPT : public MollerPlessetPT
{
public:
    UnrestrictedMollerPlessetPT(System *system, int perturbOrder, int frozenCore=0);
    virtual void solve();

private:
    UHF *m_solver;
    colvec m_fockEnergyUp;
    colvec m_fockEnergyDown;
    mat m_Cup;
    mat m_Cdown;
    int m_nElectronsUp;
    int m_nElectronsDown;

    field<mat> m_MOI_UU;
    field<mat> m_MOI_DD;
    field<mat> m_MOI_UD;

    void calc2OrderPerturb();
    void calc3OrderPerturb();
};

#endif // UNRESTRICTEDMOLLERPLESSETPT_H
