#include "rhf.h"


RHF::RHF(System* system):HartreeFock(system)
{
    m_F = zeros<mat>(m_matDim, m_matDim);
    m_C = zeros<mat>(m_matDim, m_matDim);
    m_P = zeros<mat>(m_matDim, m_matDim);
    m_fockEnergy = ones<colvec>(m_matDim)*1.0E6;
    m_nElectrons = system->getNumOfElectrons();
    if (m_nElectrons % 2 == 1){
        cout <<"Error: Cannot run Restricted Hartree-Fock on odd number of electrons." << endl;
        exit(EXIT_FAILURE);
    }
    m_restrictedFactor = 2;
}


//-------------------------------------------------------------------------------------------------------
// Solves the Hartree-Fock equations (iterated), stores the energy in double energy and the coefficients
// in vec C
void RHF::solve()
{
    // Calculate integrals
    calcIntegrals();

    // Diagonalize m_S (overlap) and calculate transformation matrix m_V
    // such that m_V.t()*S*m_V = I.
    diagOverlap();

    // Initialize density matrix (non-interacting case)
    m_P = zeros<mat>(m_matDim, m_matDim);

    // Iterate until the fock energy has converged
    double fockEnergyOld;
    double energyDiff = 1.0;
    while (energyDiff > m_toler){
        fockEnergyOld = m_fockEnergy(0);
        buildFockMatrix();
        solveSingle(m_F, m_C, m_P, m_fockEnergy, m_nElectrons);
        energyDiff = fabs(fockEnergyOld - m_fockEnergy(0));
    }
    cout << "done solve-loop" << endl;

    // Calculate energy (not equal to Fock energy)
    m_energy = 0;

    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            m_energy += 0.5*(m_h(i,j) + m_F(i,j))*m_P(i,j);
        }
    }
    m_energy += m_system->getNucleiPotential();
}

//----------------------------------------------------------------------------------------------------------------
field<mat> RHF::getCoeff(){
    field<mat> Coeff(1);
    Coeff(0) = m_C;
    return Coeff;
}

//----------------------------------------------------------------------------------------------------------------
field<mat> RHF::getDensityMatrix()
{
    field<mat> P(1);
    P(0) = m_P;
    return P;
}

//---------------------------------------------------------------------------------------------------------------
field<colvec> RHF::getFockEnergy()
{
    field<colvec> fockEnergy(1);
    fockEnergy(0) = m_fockEnergy;
    return fockEnergy;
}


//-----------------------------------------------------------------------------------------------------------------
// Builds F matrix (kinetic part, nuclear attraction part and electron-electron repulsion part, i.e.left hand side)

void RHF::buildFockMatrix()
{
    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){

            // One-electron integrals
            m_F(i,j) = m_h(i,j);

            // Add two-electron integrals
            for (int k = 0; k < m_matDim; k++){
                for (int l = 0; l < m_matDim; l++){
                    m_F(i,j) += 0.5*m_P(l,k)*(2*m_Q(i,k)(j,l) - m_Q(i,k)(l,j));
                }
            }
        }
    }
}
