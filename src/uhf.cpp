#include "uhf.h"

UHF::UHF(System *newSystem):HartreeFock(newSystem)
{
    Fup = zeros<mat>(matDim, matDim);
    Fdown = zeros<mat>(matDim, matDim);
    Cup = zeros<mat>(matDim, matDim);
    Cdown = zeros<mat>(matDim, matDim);
    Pup = zeros<mat>(matDim, matDim);
    Pup(0,1)=0.1; // Introduce an assymetry between the spin up and spin down orbitals
    Pdown = zeros<mat>(matDim, matDim);
    fockEnergyUp = ones<colvec>(matDim)*1.0E6;
    fockEnergyDown = ones<colvec>(matDim)*1.0E6;
    perturbOrder = 1;
}

UHF::UHF(System *newSystem, int newPerturbOrder):HartreeFock(newSystem)
{
    Fup = zeros<mat>(matDim, matDim);
    Fdown = zeros<mat>(matDim, matDim);
    Cup = zeros<mat>(matDim, matDim);
    Cdown = zeros<mat>(matDim, matDim);
    Pup = zeros<mat>(matDim, matDim);
    Pup(0,1)=0.1; // Introduce an assymetry between the spin up and spin down orbitals
    Pdown = zeros<mat>(matDim, matDim);
    fockEnergyUp = ones<colvec>(matDim)*1.0E6;
    fockEnergyDown = ones<colvec>(matDim)*1.0E6;
    perturbOrder = newPerturbOrder;
}





//-----------------------------------------------------------------------------------------------------------------
// Builds F+ and F- matrices (kinetic part, nuclear attraction part and electron-electron repulsion part, i.e.left hand side)

void UHF::buildFockMatrix()
{
    for (int i = 0; i < matDim; i++){
        for (int j = 0; j < matDim; j++){

            // One-electron integrals
            Fup(i,j) = h(i,j);
            Fdown(i,j) = h(i,j);

            // Add two-electron integrals
            for (int k = 0; k < matDim; k++){
                for (int l = 0; l < matDim; l++){
                    Fup(i,j) += 0.5*Pup(l,k)*(Q[i][k][j][l] - Q[i][k][l][j]) + 0.5*Pdown(l,k)*Q[i][k][j][l];
                    Fdown(i,j) += 0.5*Pdown(l,k)*(Q[i][k][j][l] - Q[i][k][l][j]) + 0.5*Pup(l,k)*Q[i][k][j][l];
                }
            }
        }
    }
}


//----------------------------------------------------------------------------------------------------------------
mat UHF::getCoeff(){
    return Cup;
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
    while (energyDiff > toler){
        fockEnergyUpOld = fockEnergyUp(0);
        fockEnergyDownOld = fockEnergyDown(0);
        buildFockMatrix();
        solveSingle(Fup, Cup, Pup, fockEnergyUp);
        solveSingle(Fdown, Cdown, Pdown, fockEnergyDown);
        energyDiff = fabs(fockEnergyUpOld + fockEnergyDownOld - fockEnergyUp(0) - fockEnergyDown(0));
    }

    // Calculate energy (not equal to Fock energy)
    energy = 0;

    for (int i = 0; i < matDim; i++){
        for (int j = 0; j < matDim; j++){
            energy += 0.25*((Pup(i,j) + Pdown(i,j))*h(i, j) + Fup(i,j)*Pup(i,j) + Fdown(i,j)*Pdown(i,j));
        }
    }
    energy += system->getNucleiPotential();

    // Perturbative terms
    if (perturbOrder == 2){
        energy += perturbation2order();
    } else if (perturbOrder == 1) {
    } else {
        cout << "Error. Only first and second order perturbation has been implemented." << endl;
        exit(EXIT_FAILURE);
    }
}


double UHF::perturbation2order(){
    int nHStates = nElectrons/2;     // later: divide between nElectronsUp and nElectronsDown
    int nPStates = matDim - nElectrons/2;
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
                    for (int p = 0; p < matDim; p++){
                        for (int q = 0; q < matDim; q++){
                            for (int r = 0; r < matDim; r++){
                                for (int s = 0; s < matDim; s++){
                                    orbitalIntegralsUU(i,j)(a,b) += Cup(p,i)*Cup(q,j)*Cup(r,a+nHStates)*Cup(s,b+nHStates)*Q[p][q][r][s];
                                    orbitalIntegralsDD(i,j)(a,b) += Cdown(p,i)*Cdown(q,j)*Cdown(r,a+nHStates)*Cdown(s,b+nHStates)*Q[p][q][r][s];
                                    orbitalIntegralsUD(i,j)(a,b) += Cup(p,i)*Cdown(q,j)*Cup(r,a+nHStates)*Cdown(s,b+nHStates)*Q[p][q][r][s];
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
                                       (fockEnergyUp(i) + fockEnergyUp(j) - fockEnergyUp(a+nHStates) - fockEnergyUp(b+nHStates));
                    energyTemp2 += 0.25*(orbitalIntegralsDD(i,j)(a,b) - orbitalIntegralsDD(i,j)(b,a))*(orbitalIntegralsDD(i,j)(a,b) - orbitalIntegralsDD(i,j)(b,a))/
                                        (fockEnergyDown(i) + fockEnergyDown(j) - fockEnergyDown(a+nHStates) - fockEnergyDown(b+nHStates));
                    energyTemp3 += orbitalIntegralsUD(i,j)(a,b)*orbitalIntegralsUD(i,j)(a,b)/
                                        (fockEnergyUp(i) + fockEnergyDown(j) - fockEnergyUp(a+nHStates) - fockEnergyDown(b+nHStates));
                }
            }
        }
    }
    energyMP2 = energyTemp1 + energyTemp2 + energyTemp3;
    return energyMP2;

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
