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
    m_perturbOrder = 1;
}

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
    m_perturbOrder = perturbOrder;
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
                    m_Fup(i,j) += 0.5*m_Pup(l,k)*(m_Q[i][k][j][l] - m_Q[i][k][l][j]) + 0.5*m_Pdown(l,k)*m_Q[i][k][j][l];
                    m_Fdown(i,j) += 0.5*m_Pdown(l,k)*(m_Q[i][k][j][l] - m_Q[i][k][l][j]) + 0.5*m_Pup(l,k)*m_Q[i][k][j][l];
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
        solveSingle(m_Fup, m_Cup, m_Pup, m_fockEnergyUp);
        solveSingle(m_Fdown, m_Cdown, m_Pdown, m_fockEnergyDown);
        energyDiff = fabs(fockEnergyUpOld + fockEnergyDownOld - m_fockEnergyUp(0) - m_fockEnergyDown(0));
    }

    // Calculate energy (not equal to Fock energy)
    m_energy = 0;

    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            m_energy += 0.25*((m_Pup(i,j) + m_Pdown(i,j))*m_h(i, j) + m_Fup(i,j)*m_Pup(i,j) + m_Fdown(i,j)*m_Pdown(i,j));
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
    int nHStates = m_nElectrons/2;     // later: divide between nElectronsUp and nElectronsDown
    int nPStates = m_matDim - m_nElectrons/2;
    field<mat> orbitalIntegralsUU(nHStates, nHStates);
    field<mat> orbitalIntegralsDD(nHStates, nHStates);
    field<mat> orbitalIntegralsUD(nHStates, nHStates);
    for (int i = 0; i < nHStates; i++){
        for (int j = 0; j < nHStates; j++){
            orbitalIntegralsUU(i,j) = zeros(nPStates, nPStates);
            orbitalIntegralsDD(i,j) = zeros(nPStates, nPStates);
            orbitalIntegralsUD(i,j) = zeros(nPStates, nPStates);
            for (int a = 0; a < nPStates; a++){
                for (int b = 0; b < nPStates; b++){
                    for (int p = 0; p < m_matDim; p++){
                        for (int q = 0; q < m_matDim; q++){
                            for (int r = 0; r < m_matDim; r++){
                                for (int s = 0; s < m_matDim; s++){
                                    orbitalIntegralsUU(i,j)(a,b) += m_Cup(p,i)*m_Cup(q,j)*m_Cup(r,a+nHStates)*m_Cup(s,b+nHStates)*m_Q[p][q][r][s];
                                    orbitalIntegralsDD(i,j)(a,b) += m_Cdown(p,i)*m_Cdown(q,j)*m_Cdown(r,a+nHStates)*m_Cdown(s,b+nHStates)*m_Q[p][q][r][s];
                                    orbitalIntegralsUD(i,j)(a,b) += m_Cup(p,i)*m_Cdown(q,j)*m_Cup(r,a+nHStates)*m_Cdown(s,b+nHStates)*m_Q[p][q][r][s];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    double energyTemp1 = 0;
    double energyTemp2 = 0;
    double energyTemp3 = 0;
    for (int i = 0; i < nHStates; i++){
        for (int j = 0; j < nHStates; j++){
            for (int a = 0; a < nPStates; a++){
                for (int b = 0; b < nPStates; b++){
                    energyTemp1 += 0.25*(orbitalIntegralsUU(i,j)(a,b) - orbitalIntegralsUU(i,j)(b,a))*(orbitalIntegralsUU(i,j)(a,b) - orbitalIntegralsUU(i,j)(b,a))/
                                       (m_fockEnergyUp(i) + m_fockEnergyUp(j) - m_fockEnergyUp(a+nHStates) - m_fockEnergyUp(b+nHStates));
                    energyTemp2 += 0.25*(orbitalIntegralsDD(i,j)(a,b) - orbitalIntegralsDD(i,j)(b,a))*(orbitalIntegralsDD(i,j)(a,b) - orbitalIntegralsDD(i,j)(b,a))/
                                        (m_fockEnergyDown(i) + m_fockEnergyDown(j) - m_fockEnergyDown(a+nHStates) - m_fockEnergyDown(b+nHStates));
                    energyTemp3 += orbitalIntegralsUD(i,j)(a,b)*orbitalIntegralsUD(i,j)(a,b)/
                                        (m_fockEnergyUp(i) + m_fockEnergyDown(j) - m_fockEnergyUp(a+nHStates) - m_fockEnergyDown(b+nHStates));
                }
            }
        }
    }
    m_energyMP2 = energyTemp1 + energyTemp2 + energyTemp3;
    return m_energyMP2;

    //    double energyTemp1, energyTemp2, energyTemp3, energyTemp4;
//    for (int i = 0; i < nElectrons/2; i++){
//        for (int j = 0; j < nElectrons/2; j++){
//            for (int a = nElectrons/2; a < matDim; a++){
//                for (int b = nElectrons/2; b < matDim; b++){
//                    energyTemp1 = 0;
//                    energyTemp2 = 0;
//                    energyTemp3 = 0;
//                    energyTemp4 = 0;
//                    if (a == b || i == j){
//                        for (int p = 0; p < matDim; p++){
//                            for (int q = 0; q < matDim; q++){
//                                for (int r = 0; r < matDim; r++){
//                                    for (int s = 0; s < matDim; s++){
//                                        energyTemp2 += Cdown(p,i)*Cdown(r,a)*Cup(q,j)*Cup(s,b)*Q[p][q][r][s];
//                                        energyTemp4 += -Cdown(p,i)*Cup(r,a)*Cup(q,j)*Cdown(s,b)*Q[p][q][s][r];
//                                    }
//                                }
//                            }
//                        }
//                    } else {
//                        for (int p = 0; p < matDim; p++){
//                            for (int q = 0; q < matDim; q++){
//                                for (int r = 0; r < matDim; r++){
//                                    for (int s = 0; s < matDim; s++){
//                                        energyTemp1 += Cup(p,i)*Cup(r,a)*Cup(q,j)*Cup(s,b)*(Q[p][q][r][s] - Q[p][q][s][r]);
//                                        energyTemp2 += Cdown(p,i)*Cdown(r,a)*Cup(q,j)*Cup(s,b)*Q[p][q][r][s];
//                                        energyTemp3 += Cdown(p,i)*Cdown(r,a)*Cdown(q,j)*Cdown(s,b)*(Q[p][q][r][s] - Q[p][q][s][r]);
//                                        energyTemp4 += -Cdown(p,i)*Cup(r,a)*Cup(q,j)*Cdown(s,b)*Q[p][q][s][r];
//                                    }
//                                }
//                            }
//                        }
//                    }
//                    energyMP2 += (energyTemp1*energyTemp1/(fockEnergyUp(i) - fockEnergyUp(a) + fockEnergyUp(j) - fockEnergyUp(b))
//                                  + 2*energyTemp2*energyTemp2/(fockEnergyDown(i) - fockEnergyDown(a) + fockEnergyUp(j) - fockEnergyUp(b))
//                                  + energyTemp3*energyTemp3/(fockEnergyDown(i) - fockEnergyDown(a) + fockEnergyDown(j) - fockEnergyDown(b))
//                                  + 2*energyTemp4*energyTemp4/(fockEnergyDown(i) - fockEnergyUp(a) + fockEnergyUp(j) - fockEnergyDown(b)))/4;
//                }
//            }
//        }
//    }
//    return energyMP2;
}
