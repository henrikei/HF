#include "rhf.h"


RHF::RHF(System* system, int perturbOrder):HartreeFock(system)
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
    m_perturbOrder = perturbOrder;
    m_energyMP2 = 0;
}


//-------------------------------------------------------------------------------------------------------
// Solves the Hartree-Fock equations (iterated), stores the energy in double energy and the coefficients
// in vec C
void RHF::solve()
{
    // Calculate integrals
    calcIntegrals();

    // Iterate until the fock energy has converged
    double fockEnergyOld;
    double energyDiff = 1.0;
    int counter = 0;
    while (energyDiff > m_toler){
        counter += 1;
        fockEnergyOld = m_fockEnergy(0);
        buildFockMatrix();
        solveSingle(m_F, m_C, m_P, m_fockEnergy, m_nElectrons);
        energyDiff = fabs(fockEnergyOld - m_fockEnergy(0));
    }
    cout << counter << endl;

    // Calculate energy (not equal to Fock energy)
    m_energy = 0;

    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            m_energy += 0.5*(m_h(i,j) + m_F(i,j))*m_P(i,j);
        }
    }
    m_energy += m_system->getNucleiPotential();

    // Perturbative terms
    if (m_perturbOrder == 2){
        m_energy += perturbation2order();
    } else if (m_perturbOrder == 1) {
    } else {
        cout << "Error. Only first and second order perturbation has been implemented." << endl;
        exit(EXIT_FAILURE);
    }
}

//----------------------------------------------------------------------------------------------------------------
mat RHF::getCoeff(){
    return m_C;
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


//--------------------------------------------------------------------------------------------------------------------
// Second order perturbation
double RHF::perturbation2order(){
    int nHStates = m_nElectrons/2;
    int nPStates = m_matDim - m_nElectrons/2;
    // temp1, temp2, temp3 are temporary fields used to get from
    // Atomic Orbital Integrals to Molecular Orbital Integrals
    field<mat> temp1(nHStates, m_matDim);
    field<mat> temp2(nHStates, nHStates);
    field<mat> temp3(nHStates, nHStates);
    // orbitalIntegral = <ij|g|ab>
    field<mat> MOI(nHStates, nHStates);
    for (int i = 0; i < nHStates; i++){
        for (int q = 0; q < m_matDim; q++){
            temp1(i,q) = zeros(m_matDim, m_matDim);
        }
    }
    for (int i = 0; i < nHStates; i++){
        for (int j = 0; j < nHStates; j++){
            temp2(i,j) = zeros(m_matDim, m_matDim);
            temp3(i,j) = zeros(nPStates, m_matDim);
            MOI(i,j) = zeros(nPStates, nPStates);
        }
    }

    // Transform from Atomic Orbital Integrals to Molecular Orbital Integrals
    AOItoMOI(temp1, m_Q, m_C, 0);
    AOItoMOI(temp2, temp1, m_C, 1);
    AOItoMOI(temp3, temp2, m_C.cols(nHStates, m_matDim-1), 2);
    AOItoMOI(MOI, temp3, m_C.cols(nHStates, m_matDim-1), 3);

    // Sum up energy tems
    for (int i = 0; i < nHStates; i++){
        for (int j = 0; j < nHStates; j++){
            for (int a = 0; a < nPStates; a++){
                for (int b = 0; b < nPStates; b++){
                    m_energyMP2 += MOI(i,j)(a,b)*(2*MOI(i,j)(a,b) - MOI(j,i)(a,b))
                                /(m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(a+nHStates) - m_fockEnergy(b+nHStates));
                }
            }
        }
    }

    return m_energyMP2;
}
