#ifndef UMF_H
#define UMP_H

#include <armadillo>
#include "mollerplesset.h"
#include "hartreefock/uhf.h"

using namespace std;
using namespace arma;


class UMP : public MollerPlesset
{
public:
    UMP(System *system, int perturbOrder, int frozenCore=0);
    virtual void solve();
    virtual field<mat> getDensityMatrix();

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

#endif // UMP_H
