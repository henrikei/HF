#include "uhf.h"

UHF::UHF(System *system, int perturbOrder):HartreeFock(system)
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
    m_perturbOrder = perturbOrder;

    m_MOI_UU.set_size(m_matDim, m_matDim);
    m_MOI_DD.set_size(m_matDim, m_matDim);
    m_MOI_UD.set_size(m_matDim, m_matDim);
    for(int i = 0; i < m_matDim; i++){
        for(int j = 0; j < m_matDim; j++){
            m_MOI_UU(i,j) = zeros(m_matDim, m_matDim);
            m_MOI_DD(i,j) = zeros(m_matDim, m_matDim);
            m_MOI_UD(i,j) = zeros(m_matDim, m_matDim);
        }
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
mat UHF::getCoeff(){
    return m_Cup;
}


//-------------------------------------------------------------------------------------------------------
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

    // Perturbative terms
    if (m_perturbOrder == 2){
        m_energy += perturbation2order();
    } else if (m_perturbOrder == 1) {
    } else {
        cout << "Error. Only first and second order perturbation has been implemented." << endl;
        exit(EXIT_FAILURE);
    }
}


double UHF::perturbation2order(){
    field<mat> tempUU1(m_matDim, m_matDim);
    field<mat> tempDD1(m_matDim, m_matDim);
    field<mat> tempUD1(m_matDim, m_matDim);
    field<mat> tempUU2(m_matDim, m_matDim);
    field<mat> tempDD2(m_matDim, m_matDim);
    field<mat> tempUD2(m_matDim, m_matDim);
    field<mat> tempUU3(m_matDim, m_matDim);
    field<mat> tempDD3(m_matDim, m_matDim);
    field<mat> tempUD3(m_matDim, m_matDim);
    // MOI = Molecular Orbital Integral
    // AOI = Atomic Orbital Integral

    for(int i = 0; i < m_matDim; i++){
        for(int j = 0; j < m_matDim; j++){
            tempUU1(i,j) = zeros(m_matDim, m_matDim);
            tempDD1(i,j) = zeros(m_matDim, m_matDim);
            tempUD1(i,j) = zeros(m_matDim, m_matDim);
            tempUU2(i,j) = zeros(m_matDim, m_matDim);
            tempDD2(i,j) = zeros(m_matDim, m_matDim);
            tempUD2(i,j) = zeros(m_matDim, m_matDim);
            tempUU3(i,j) = zeros(m_matDim, m_matDim);
            tempDD3(i,j) = zeros(m_matDim, m_matDim);
            tempUD3(i,j) = zeros(m_matDim, m_matDim);
        }
    }

    // Up-Up AOI to Up-Up MOI
    AOItoMOI(tempUU1, m_Q, m_Cup, 0);
    AOItoMOI(tempUU2, tempUU1, m_Cup, 1);
    AOItoMOI(tempUU3, tempUU2, m_Cup, 2);
    AOItoMOI(m_MOI_UU, tempUU3, m_Cup, 3);

    // Down-Down AOI to Down-Down MOI
    AOItoMOI(tempDD1, m_Q, m_Cdown, 0);
    AOItoMOI(tempDD2, tempDD1, m_Cdown, 1);
    AOItoMOI(tempDD3, tempDD2, m_Cdown, 2);
    AOItoMOI(m_MOI_DD, tempDD3, m_Cdown, 3);

    // Up-Down AOI to Up-Down MOI
    AOItoMOI(tempUD1, m_Q, m_Cup, 0);
    AOItoMOI(tempUD2, tempUD1, m_Cdown, 1);
    AOItoMOI(tempUD3, tempUD2, m_Cup, 2);
    AOItoMOI(m_MOI_UD, tempUD3, m_Cdown, 3);


    // Sum up energy terms
    for (int i = 0; i < m_nElectronsUp; i++){
        for (int j = 0; j < m_nElectronsUp; j++){
            for (int a = m_nElectronsUp; a < m_matDim; a++){
                for (int b = m_nElectronsUp; b < m_matDim; b++){
                    m_energyMP2 += 0.25*(m_MOI_UU(i,j)(a,b) - m_MOI_UU(i,j)(b,a))*(m_MOI_UU(i,j)(a,b) - m_MOI_UU(i,j)(b,a))/
                            (m_fockEnergyUp(i) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyUp(b));
                }
            }
        }
    }
    for (int i = 0; i < m_nElectronsDown; i++){
        for (int j = 0; j < m_nElectronsDown; j++){
            for (int a = m_nElectronsDown; a < m_matDim; a++){
                for (int b = m_nElectronsDown; b < m_matDim; b++){
                    m_energyMP2 += 0.25*(m_MOI_DD(i,j)(a,b) - m_MOI_DD(i,j)(b,a))*(m_MOI_DD(i,j)(a,b) - m_MOI_DD(i,j)(b,a))/
                            (m_fockEnergyDown(i) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyDown(b));
                }
            }
        }
    }
    for (int i = 0; i < m_nElectronsUp; i++){
        for (int j = 0; j < m_nElectronsDown; j++){
            for (int a = m_nElectronsUp; a < m_matDim; a++){
                for (int b = m_nElectronsDown; b < m_matDim; b++){
                    m_energyMP2 += m_MOI_UD(i,j)(a,b)*m_MOI_UD(i,j)(a,b)/
                            (m_fockEnergyUp(i) + m_fockEnergyDown(j) - m_fockEnergyUp(a) - m_fockEnergyDown(b));
                }
            }
        }
    }
    return m_energyMP2;
}

double UHF::perturbation3order()
{
    return 0;
}
