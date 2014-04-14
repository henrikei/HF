#ifndef MOLLERPLESSET_H
#define MOLLERPLESSET_H

#include <armadillo>
#include <mpi.h>
#include <system/system.h>
#include <hartreefock/hartreefock.h>

#define RUN_MPI

class MollerPlesset
{
public:
    MollerPlesset(System *system, int perturbOrder, int frozenCore);
    virtual void solve()=0;
    double getEnergyHF();
    double getEnergy2order();
    double getEnergy3order();

protected:
    System *m_system;
    field<mat> m_AOI;

    int m_perturbOrder;
    int m_matDim;
    double m_energyHF;
    double m_energy2Order;
    double m_energy3Order;

    int m_frozenCore;

    void AOItoMOI(field<mat> &MOI, field<mat> AOI, mat C, int index);
};

#endif // MOLLERPLESSET_H
