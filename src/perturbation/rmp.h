#ifndef RMP_H
#define RMP_H

#include <armadillo>
#include <time.h>
#include <system/system.h>
#include <perturbation/mollerplesset.h>
#include <hartreefock/rhf.h>

using namespace std;
using namespace arma;

class RMP : public MollerPlesset
{
public:
    RMP(System *system, int perturbOrder, int frozenCore=0);
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

#endif // RMP_H
