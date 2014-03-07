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
    m_restrictedFactor = 1;
    m_nElectronsUp = system->getNumOfElectrons()/2;
    m_nElectronsDown = m_nElectronsUp;
    if (system->getNumOfElectrons() % 2 == 1){
        m_nElectronsUp += 1;
    }
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
    int nHStatesUp = m_nElectronsUp;
    int nHStatesDown = m_nElectronsDown;
    int nPStatesUp = m_matDim - m_nElectronsUp;
    int nPStatesDown = m_matDim - m_nElectronsDown;
    field<mat> tempUU1(nHStatesUp, m_matDim);
    field<mat> tempDD1(nHStatesDown, m_matDim);
    field<mat> tempUD1(nHStatesUp, m_matDim);
    field<mat> tempUU2(nHStatesUp, nHStatesUp);
    field<mat> tempDD2(nHStatesDown, nHStatesDown);
    field<mat> tempUD2(nHStatesUp, nHStatesDown);
    field<mat> tempUU3(nHStatesUp, nHStatesUp);
    field<mat> tempDD3(nHStatesDown, nHStatesDown);
    field<mat> tempUD3(nHStatesUp, nHStatesDown);
    // MOI = Molecular Orbital Integral
    // AOI = Atomic Orbital Integral
    field<mat> MOI_UU(nHStatesUp, nHStatesUp);
    field<mat> MOI_DD(nHStatesDown, nHStatesDown);
    field<mat> MOI_UD(nHStatesUp, nHStatesDown);

    // Up-Up AOI to Up-Up MOI
    for(int i = 0; i < nHStatesUp; i++){
        for(int j = 0; j < m_matDim; j++){
            tempUU1(i,j) = zeros(m_matDim, m_matDim);
        }
    }
    for(int i = 0; i < nHStatesUp; i++){
        for(int j = 0; j < nHStatesUp; j++){
            tempUU2(i,j) = zeros(m_matDim, m_matDim);
            tempUU3(i,j) = zeros(nPStatesUp, m_matDim);
            MOI_UU(i,j) = zeros(nPStatesUp, nPStatesUp);
        }
    }
    AOItoMOI(tempUU1, m_Q, m_Cup, 0);
    AOItoMOI(tempUU2, tempUU1, m_Cup, 1);
    AOItoMOI(tempUU3, tempUU2, m_Cup.cols(nHStatesUp, m_matDim-1), 2);
    AOItoMOI(MOI_UU, tempUU3, m_Cup.cols(nHStatesUp, m_matDim-1), 3);

    // Down-Down AOI to Down-Down MOI
    for(int i = 0; i < nHStatesDown; i++){
        for(int j = 0; j < m_matDim; j++){
            tempDD1(i,j) = zeros(m_matDim, m_matDim);
        }
    }
    for(int i = 0; i < nHStatesDown; i++){
        for(int j = 0; j < nHStatesDown; j++){
            tempDD2(i,j) = zeros(m_matDim, m_matDim);
            tempDD3(i,j) = zeros(nPStatesDown, m_matDim);
            MOI_DD(i,j) = zeros(nPStatesDown, nPStatesDown);
        }
    }
    AOItoMOI(tempDD1, m_Q, m_Cdown, 0);
    AOItoMOI(tempDD2, tempDD1, m_Cdown, 1);
    AOItoMOI(tempDD3, tempDD2, m_Cdown.cols(nHStatesDown, m_matDim-1), 2);
    AOItoMOI(MOI_DD, tempDD3, m_Cdown.cols(nHStatesDown, m_matDim-1), 3);

    // Up-Down AOI to Up-Down MOI
    for(int i = 0; i < nHStatesUp; i++){
        for(int j = 0; j < m_matDim; j++){
            tempUD1(i,j) = zeros(m_matDim, m_matDim);
        }
    }
    for(int i = 0; i < nHStatesUp; i++){
        for(int j = 0; j < nHStatesDown; j++){
            tempUD2(i,j) = zeros(m_matDim, m_matDim);
            tempUD3(i,j) = zeros(nPStatesUp, m_matDim);
            MOI_UD(i,j) = zeros(nPStatesUp, nPStatesDown);
        }
    }
    AOItoMOI(tempUD1, m_Q, m_Cup, 0);
    AOItoMOI(tempUD2, tempUD1, m_Cdown, 1);
    AOItoMOI(tempUD3, tempUD2, m_Cup.cols(nHStatesUp, m_matDim-1), 2);
    AOItoMOI(MOI_UD, tempUD3, m_Cdown.cols(nHStatesDown, m_matDim-1), 3);


    // Sum up energy terms
    for (int i = 0; i < nHStatesUp; i++){
        for (int j = 0; j < nHStatesUp; j++){
            for (int a = 0; a < nPStatesUp; a++){
                for (int b = 0; b < nPStatesUp; b++){
                    m_energyMP2 += 0.25*(MOI_UU(i,j)(a,b) - MOI_UU(i,j)(b,a))*(MOI_UU(i,j)(a,b) - MOI_UU(i,j)(b,a))/
                            (m_fockEnergyUp(i) + m_fockEnergyUp(j) - m_fockEnergyUp(a+nHStatesUp) - m_fockEnergyUp(b+nHStatesUp));
                }
            }
        }
    }
    for (int i = 0; i < nHStatesDown; i++){
        for (int j = 0; j < nHStatesDown; j++){
            for (int a = 0; a < nPStatesDown; a++){
                for (int b = 0; b < nPStatesDown; b++){
                    m_energyMP2 += 0.25*(MOI_DD(i,j)(a,b) - MOI_DD(i,j)(b,a))*(MOI_DD(i,j)(a,b) - MOI_DD(i,j)(b,a))/
                            (m_fockEnergyDown(i) + m_fockEnergyDown(j) - m_fockEnergyDown(a+nHStatesDown) - m_fockEnergyDown(b+nHStatesDown));
                }
            }
        }
    }
    for (int i = 0; i < nHStatesUp; i++){
        for (int j = 0; j < nHStatesDown; j++){
            for (int a = 0; a < nPStatesUp; a++){
                for (int b = 0; b < nPStatesDown; b++){
                    m_energyMP2 += MOI_UD(i,j)(a,b)*MOI_UD(i,j)(a,b)/
                            (m_fockEnergyUp(i) + m_fockEnergyDown(j) - m_fockEnergyUp(a+nHStatesUp) - m_fockEnergyDown(b+nHStatesDown));
                }
            }
        }
    }
    return m_energyMP2;
}

void UHF::AOItoMOI(field<mat>& MOI, field<mat> AOI, mat C, int index)
{
    int I = MOI.n_rows;
    int J = MOI.n_cols;
    int K = MOI(0,0).n_rows;
    int L = MOI(0,0).n_cols;

    if (index == 0){
        for (int i = 0; i < I; i++){
            for (int j = 0; j < J; j++){
                for (int k = 0; k < K; k++){
                    for (int l = 0; l < L; l++){
                        for (int m = 0; m < m_matDim; m++){
                            MOI(i,j)(k,l) += C(m,i)*AOI(m,j)(k,l);
                        }
                    }
                }
            }
        }
    } else if (index == 1){
        for (int i = 0; i < I; i++){
            for (int j = 0; j < J; j++){
                for (int k = 0; k < K; k++){
                    for (int l = 0; l < L; l++){
                        for (int m = 0; m < m_matDim; m++){
                            MOI(i,j)(k,l) += C(m,j)*AOI(i,m)(k,l);
                        }
                    }
                }
            }
        }
    } else if (index == 2){
        for (int i = 0; i < I; i++){
            for (int j = 0; j < J; j++){
                for (int k = 0; k < K; k++){
                    for (int l = 0; l < L; l++){
                        for (int m = 0; m < m_matDim; m++){
                            MOI(i,j)(k,l) += C(m,k)*AOI(i,j)(m,l);
                        }
                    }
                }
            }
        }
    } else if (index == 3){
        for (int i = 0; i < I; i++){
            for (int j = 0; j < J; j++){
                for (int k = 0; k < K; k++){
                    for (int l = 0; l < L; l++){
                        for (int m = 0; m < m_matDim; m++){
                            MOI(i,j)(k,l) += C(m,l)*AOI(i,j)(k,m);
                        }
                    }
                }
            }
        }
    }

    //    double liste[4];
    //    double listeC[4];
    //    for (int i = 0; i < I; i++){
    //        for (int j = 0; j < J; j++){
    //            for (int k = 0; k < K; k++){
    //                for (int l = 0; l < L; l++){
    //                    for (int m = 0; m < m_matDim; m++){
    //                        liste[0] = AOI(m,j)(k,l);
    //                        liste[1] = AOI(i,m)(k,l);
    //                        liste[2] = AOI(i,j)(m,l);
    //                        liste[3] = AOI(i,j)(k,m);
    //                        listeC[0] = C(i,m);
    //                        listeC[1] = C(j,m);
    //                        listeC[2] = C(k,m);
    //                        listeC[3] = C(l,m);
    //                        MOI(i,j)(k,l) += listeC[index] * liste[index];
    //                    }
    //                }
    //            }
    //        }
    //    }
}
