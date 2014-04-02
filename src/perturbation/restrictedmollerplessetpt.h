#ifndef RESTRICTEDMOLLERPLESSETPT_H
#define RESTRICTEDMOLLERPLESSETPT_H

#include <armadillo>
#include <system.h>
#include <perturbation/mollerplessetpt.h>
#include <rhf.h>

using namespace std;
using namespace arma;

class RestrictedMollerPlessetPT : public MollerPlessetPT
{
public:
    RestrictedMollerPlessetPT(System *system, int perturbOrder, int frozenCore=0);
    virtual void solve();

private:
    RHF *m_solver;
    colvec m_fockEnergy;
    mat m_C;
    int m_nElectrons;

    field<mat> m_MOI;

    void calc2OrderPerturb();
    void calc3OrderPerturb();
};

#endif // RESTRICTEDMOLLERPLESSETPT_H
