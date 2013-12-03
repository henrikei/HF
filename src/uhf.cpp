#include "uhf.h"

UHF::UHF(System *newSystem):HartreeFock(newSystem)
{
    Fup = zeros<mat>(matDim, matDim);
    Fdown = zeros<mat>(matDim, matDim);
    Cup = zeros<mat>(matDim, nElectrons/2);
    Cdown = zeros<mat>(matDim, nElectrons/2);
    Pup = zeros<mat>(matDim, matDim);
    Pdown = zeros<mat>(matDim, matDim);
    fockEnergyUp = 1.0E6;
    fockEnergyDown = 1.0E6;
}




//-----------------------------------------------------------------------------------------------------------------
// Builds F+ and F- matrices (kinetic part, nuclear attraction part and electron-electron repulsion part, i.e.left hand side)

void UHF::buildMatrix()
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
    Cup = zeros<mat>(matDim, nElectrons/2);
    Cdown = zeros<mat>(matDim, nElectrons/2);

    calcIntegrals();

    // Iterate until the fock energy has converged
    while (energyDiff > toler){
        fockEnergyUpOld = fockEnergyUp;
        fockEnergyDownOld = fockEnergyDown;
        buildMatrix();
        solveSingle(Fup, Cup, Pup, fockEnergyUp);
        solveSingle(Fdown, Cdown, Pdown, fockEnergyDown);
        energyDiff = fabs(fockEnergyUpOld + fockEnergyDownOld - fockEnergyUp - fockEnergyDown);
    }

    // Calculate energy (not equal to Fock energy)
    energy = 0;

    for (int i = 0; i < matDim; i++){
        for (int j = 0; j < matDim; j++){
            energy += 0.5*(Pup(i,j) + Pdown(i,j))*h(i, j);
            for (int k = 0; k < matDim; k++){
                for (int l = 0; l < matDim; l++){
                    energy += 0.25*Pup(i,j)*Pup(l,k)*(Q[i][k][j][l] - Q[i][k][l][j]) + 0.25*Pup(i,j)*Pdown(l,k)*Q[i][k][j][l]
                            + 0.25*Pdown(i,j)*Pdown(l,k)*(Q[i][k][j][l] - Q[i][k][l][j]) + 0.25*Pdown(i,j)*Pup(l,k)*Q[i][k][j][l];
                }
            }
        }
    }
    energy += system->getNucleiPotential();
}
