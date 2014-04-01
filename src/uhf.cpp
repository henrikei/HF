#include "uhf.h"

UHF::UHF(System *system):HartreeFock(system)
{
    m_Fup = zeros<mat>(m_matDim, m_matDim);
    m_Fdown = zeros<mat>(m_matDim, m_matDim);
    m_Cup = zeros<mat>(m_matDim, m_matDim);
    m_Cdown = zeros<mat>(m_matDim, m_matDim);
    m_Pup = zeros<mat>(m_matDim, m_matDim);
    m_Pup(0,1)=0.1; // Introduce an assymetry between the spin up and spin down orbitals
    m_Pdown = zeros<mat>(m_matDim, m_matDim);
    m_fockEnergyUp = ones<colvec>(m_matDim)*1.0E6;
    m_fockEnergyDown = ones<colvec>(m_matDim)*1.0E6;
    m_restrictedFactor = 1;
    m_nElectronsUp = system->getNumOfElectrons()/2;
    m_nElectronsDown = m_nElectronsUp;
    if (system->getNumOfElectrons() % 2 == 1){
        m_nElectronsUp += 1;
    }
}



//-----------------------------------------------------------------------------------------------------------------
// Builds F+ and F- matrices (kinetic part, nuclear attraction part and electron-electron repulsion part, i.e.left hand side)

void UHF::buildFockMatrix()
{
    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){

            // One-electron integrals
            m_Fup(i,j) = m_h(i,j);
            m_Fdown(i,j) = m_h(i,j);

            // Add two-electron integrals
            for (int k = 0; k < m_matDim; k++){
                for (int l = 0; l < m_matDim; l++){
                    m_Fup(i,j) += m_Pup(l,k)*(m_Q(i,k)(j,l) - m_Q(i,k)(l,j)) + m_Pdown(l,k)*m_Q(i,k)(j,l);
                    m_Fdown(i,j) += m_Pdown(l,k)*(m_Q(i,k)(j,l) - m_Q(i,k)(l,j)) + m_Pup(l,k)*m_Q(i,k)(j,l);
                }
            }
        }
    }
}


//----------------------------------------------------------------------------------------------------------------
field<mat> UHF::getCoeff(){
    field<mat> Coeff(2);
    Coeff(0) = m_Cup;
    Coeff(1) = m_Cdown;
    return Coeff;
}

//------------------------------------------------------------------------------------------------------------------
field<colvec> UHF::getFockEnergy()
{
    field<colvec> fockEnergy(2);
    fockEnergy(0) = m_fockEnergyUp;
    fockEnergy(1) = m_fockEnergyDown;
    return fockEnergy;
}

int UHF::getNumOfElectronsUp()
{
    return m_nElectronsUp;
}

int UHF::getNumOfElecrtonsDown()
{
    return m_nElectronsDown;
}


//-----------------------------------------------------------------------------------------------------------------
// Solves the Hartree-Fock equations (iterated), stores the energy in double energy and the coefficients
// in vec C
void UHF::solve()
{
    double fockEnergyUpOld;
    double fockEnergyDownOld;
    double energyDiff = 1.0;

    // Calculate integrals
    calcIntegrals();

    // Iterate until the fock energy has converged
    while (energyDiff > m_toler){
        fockEnergyUpOld = m_fockEnergyUp(0);
        fockEnergyDownOld = m_fockEnergyDown(0);
        buildFockMatrix();
        solveSingle(m_Fup, m_Cup, m_Pup, m_fockEnergyUp, m_nElectronsUp);
        solveSingle(m_Fdown, m_Cdown, m_Pdown, m_fockEnergyDown, m_nElectronsDown);
        energyDiff = fabs(fockEnergyUpOld - m_fockEnergyUp(0))+ fabs(fockEnergyDownOld - m_fockEnergyDown(0));
    }

    // Calculate energy (not equal to Fock energy)
    m_energy = 0;

    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            m_energy += 0.5*((m_Pup(i,j) + m_Pdown(i,j))*m_h(i, j) + m_Fup(i,j)*m_Pup(i,j) + m_Fdown(i,j)*m_Pdown(i,j));
        }
    }
    m_energy += m_system->getNucleiPotential();
}
