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
    m_nElectronsUp = system->getNumOfElectronsUp();
    m_nElectronsDown = system->getNumOfElectronsDown();
}



//-----------------------------------------------------------------------------------------------------------------
// Builds F+ and F- matrices (kinetic part, nuclear attraction part and electron-electron repulsion part, i.e.left hand side)

void UHF::buildFockMatrix()
{
    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < i+1; j++){

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
            m_Fup(j,i) = m_Fup(i,j);
            m_Fdown(j,i) = m_Fdown(i,j);
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

field<mat> UHF::getDensityMatrix()
{
    field<mat> P(2);
    P(0) = m_Pup;
    P(1) = m_Pdown;
    return P;
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

double UHF::getSpinExpectation()
{
    double value, value_temp;
    value = 0;
    for (int i = 0; i < m_nElectronsUp; i++){
        for (int j = 0; j < m_nElectronsDown; j++){
            value_temp = 0;
            for (uint mu = 0; mu < m_Cup.n_rows; mu++){
                for (uint nu = 0; nu < m_Cdown.n_rows; nu++){
                    value_temp += m_Cup(mu,i)*m_Cdown(nu,j)*m_S(mu,nu);
                }
            }
            value -= value_temp*value_temp;
        }
    }
    double factor = (double)(m_nElectronsUp - m_nElectronsDown);
    factor /= 2;
    value += factor*(factor + 1) + (double)m_nElectronsDown;
    return value;
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

    // Diagonalize m_S (overlap) and calculate transformation matrix m_V
    // such that m_V.t()*S*m_V = I.
    diagOverlap();

    // Initialize density matrices (non-interacting initial condition)
    m_Pup = zeros<mat>(m_matDim, m_matDim); m_Pup(0,1) = 0.1;
    m_Pdown = zeros<mat>(m_matDim, m_matDim);

    // Iterate until the fock energy has converged
    while (energyDiff > m_toler){
        fockEnergyUpOld = m_fockEnergyUp(0);
        fockEnergyDownOld = m_fockEnergyDown(0);
        buildFockMatrix();
        solveSingle(m_Fup, m_Cup, m_Pup, m_fockEnergyUp, m_nElectronsUp);
        solveSingle(m_Fdown, m_Cdown, m_Pdown, m_fockEnergyDown, m_nElectronsDown);
        energyDiff = fabs(fockEnergyUpOld - m_fockEnergyUp(0))+ fabs(fockEnergyDownOld - m_fockEnergyDown(0));
    }

    cout << "done solve-loop" << endl;

    // Calculate energy (not equal to Fock energy)
    m_energy = 0;

    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            m_energy += 0.5*((m_Pup(i,j) + m_Pdown(i,j))*m_h(i, j) + m_Fup(i,j)*m_Pup(i,j) + m_Fdown(i,j)*m_Pdown(i,j));
        }
    }
    m_energy += m_system->getNucleiPotential();
}
