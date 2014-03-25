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

    m_MOI.set_size(m_matDim, m_matDim);
    for(int i = 0; i < m_matDim; i++){
        for(int j = 0; j < m_matDim; j++){
            m_MOI(i,j) = zeros(m_matDim, m_matDim);
        }
    }
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
    while (energyDiff > m_toler){
        fockEnergyOld = m_fockEnergy(0);
        buildFockMatrix();
        solveSingle(m_F, m_C, m_P, m_fockEnergy, m_nElectrons);
        energyDiff = fabs(fockEnergyOld - m_fockEnergy(0));
    }

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
    } else if (m_perturbOrder ==3) {
        m_energy += perturbation2order();
        m_energy += perturbation3order();
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
    // temp1, temp2, temp3 are temporary fields used to get from
    // Atomic Orbital Integrals to Molecular Orbital Integrals
    field<mat> temp1(m_matDim, m_matDim);
    field<mat> temp2(m_matDim, m_matDim);
    field<mat> temp3(m_matDim, m_matDim);
    // orbitalIntegral = <ij|g|ab>
    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            temp1(i,j) = zeros(m_matDim, m_matDim);
            temp2(i,j) = zeros(m_matDim, m_matDim);
            temp3(i,j) = zeros(m_matDim, m_matDim);
        }
    }

    // Transform from Atomic Orbital Integrals to Molecular Orbital Integrals
    AOItoMOI(temp1, m_Q, m_C, 0);
    AOItoMOI(temp2, temp1, m_C, 1);
    AOItoMOI(temp3, temp2, m_C, 2);
    AOItoMOI(m_MOI, temp3, m_C, 3);

    // Sum up energy tems
    for (int i = 0; i < m_nElectrons/2; i++){
        for (int j = 0; j < m_nElectrons/2; j++){
            for (int a = m_nElectrons/2; a < m_matDim; a++){
                for (int b = m_nElectrons/2; b < m_matDim; b++){
                    m_energyMP2 += m_MOI(i,j)(a,b)*(2*m_MOI(i,j)(a,b) - m_MOI(j,i)(a,b))
                                /(m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(a) - m_fockEnergy(b));
                }
            }
        }
    }

    return m_energyMP2;
}

double RHF::perturbation3order()
{
    // Contribution from particle ladder diagram
    for (int i = 0; i < m_nElectrons/2; i++){
        for (int j = 0; j < m_nElectrons/2; j++){
            for (int a = m_nElectrons/2; a < m_matDim; a++){
                for (int b = m_nElectrons/2; b < m_matDim; b++){
                    for (int c = m_nElectrons/2; c < m_matDim; c++){
                        for (int d = m_nElectrons/2; d < m_matDim; d++){
                            m_energyMP3 += (m_MOI(i,j)(a,b)*m_MOI(a,b)(c,d)*(2*m_MOI(c,d)(i,j) + 2*m_MOI(d,c)(j,i) - m_MOI(c,d)(j,i) - m_MOI(d,c)(i,j))
                                         + m_MOI(i,j)(a,b)*m_MOI(b,a)(c,d)*(2*m_MOI(c,d)(j,i) + 2*m_MOI(d,c)(i,j) - m_MOI(c,d)(i,j) - m_MOI(d,c)(j,i)))
                                          /(4*(m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(a) - m_fockEnergy(b))
                                             *(m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(c) - m_fockEnergy(d)));
                        }
                    }
                }
            }
        }
    }

    // Contribution from hole ladder diagram
    for (int i = 0; i < m_nElectrons/2; i++){
        for (int j = 0; j < m_nElectrons/2; j++){
            for (int k = 0; k < m_nElectrons/2; k++){
                for (int l = 0; l < m_nElectrons/2; l++){
                    for (int a = m_nElectrons/2; a < m_matDim; a++){
                        for (int b = m_nElectrons/2; b < m_matDim; b++){
                            m_energyMP3 += (m_MOI(i,j)(a,b)*m_MOI(a,b)(k,l)*(2*m_MOI(k,l)(i,j) + 2*m_MOI(l,k)(j,i) - m_MOI(k,l)(j,i) - m_MOI(l,k)(i,j))
                                         + m_MOI(i,j)(a,b)*m_MOI(b,a)(k,l)*(2*m_MOI(k,l)(j,i) + 2*m_MOI(l,k)(i,j) - m_MOI(k,l)(i,j) - m_MOI(l,k)(j,i)))
                                           /(4*(m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(a) - m_fockEnergy(b))
                                              *(m_fockEnergy(k) + m_fockEnergy(l) - m_fockEnergy(a) - m_fockEnergy(b)));
                        }
                    }
                }
            }
        }
    }


    // Contribution from ring diagram
    for (int i = 0; i < m_nElectrons/2; i++){
        for (int j = 0; j < m_nElectrons/2; j++){
            for (int k = 0; k < m_nElectrons/2; k++){
                for (int a = m_nElectrons/2; a < m_matDim; a++){
                    for (int b = m_nElectrons/2; b < m_matDim; b++){
                        for (int c = m_nElectrons/2; c < m_matDim; c++){
                            m_energyMP3 += -2*(m_MOI(i,j)(a,b)*m_MOI(k,b)(i,c)*(2*m_MOI(a,c)(k,j) - m_MOI(a,c)(j,k))
                                              +m_MOI(i,j)(a,b)*m_MOI(k,b)(c,i)*(2*m_MOI(a,c)(j,k) - m_MOI(a,c)(k,j))
                                              +m_MOI(i,j)(b,a)*m_MOI(k,b)(i,c)*(2*m_MOI(a,c)(j,k) - m_MOI(a,c)(k,j))
                                              +m_MOI(i,j)(b,a)*m_MOI(k,b)(c,i)*(2*m_MOI(a,c)(k,j) - 4*m_MOI(a,c)(j,k)))
                                            /((m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(a) - m_fockEnergy(b))
                                             *(m_fockEnergy(k) + m_fockEnergy(j) - m_fockEnergy(a) - m_fockEnergy(c)));
                        }
                    }
                }
            }
        }
    }

    return m_energyMP3;
}
