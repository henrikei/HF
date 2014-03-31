#include "mollerplessetpt.h"

MollerPlessetPT::MollerPlessetPT(System *system, int perturbOrder)
{
    m_system = system;
    m_AOI.set_size(1,1);
    m_AOI(0,0) = zeros(1,1);
    m_perturbOrder = perturbOrder;
    m_matDim = m_system->getTotalNumOfBasisFunc();
    m_energyHF = 0;
    m_energy2Order = 0;
    m_energy3Order = 0;
}

double MollerPlessetPT::getEnergyHF()
{
    return m_energyHF;
}

double MollerPlessetPT::getEnergy2order()
{
    return m_energy2Order;
}

double MollerPlessetPT::getEnergy3order()
{
    return m_energy3Order;
}

// Transforms Atomic Orbital Integrals to Molecular Orbital Integrals one index at a time
void MollerPlessetPT::AOItoMOI(field<mat>& MOI, field<mat> AOI, mat C, int index)
{
    int a, b, c, d, e;

    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            for (int k = 0; k < m_matDim; k++){
                for (int l = 0; l < m_matDim; l++){
                    for (int m = 0; m < m_matDim; m++){
                        if (index == 0)       { e = i; a = m; b = j; c = k; d = l;}
                        else if(index == 1)   { e = j; a = i; b = m; c = k; d = l;}
                        else if(index == 2)   { e = k; a = i; b = j; c = m; d = l;}
                        else if(index == 3)   { e = l; a = i; b = j; c = k; d = m;}
                        MOI(i,j)(k,l) += C(m,e)*AOI(a,b)(c,d);
                    }
                }
            }
        }
    }
}
