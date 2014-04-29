#ifndef MOLLERPLESSET_H
#define MOLLERPLESSET_H

#include <armadillo>
#include <system/system.h>
#include <hartreefock/hartreefock.h>
#ifdef RUN_MPI
#include <mpi.h>
#endif


class MollerPlesset
{
public:
    MollerPlesset(System *system, int perturbOrder, int frozenCore);
    virtual void solve()=0;
    double getEnergy();
    double getEnergyHF();
    double getEnergy2order();
    double getEnergy3order();
    virtual field<mat> getDensityMatrix()=0;

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
