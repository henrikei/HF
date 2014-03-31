#ifndef MOLLERPLESSETPT_H
#define MOLLERPLESSETPT_H

#include <armadillo>
#include <system.h>
#include <hartreefock.h>

class MollerPlessetPT
{
public:
    MollerPlessetPT(System *system, int perturbOrder);
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

    void AOItoMOI(field<mat> &MOI, field<mat> AOI, mat C, int index);
};

#endif // MOLLERPLESSETPT_H
