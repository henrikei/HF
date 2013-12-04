#include "rhf.h"

RHF::RHF(System* newSystem):HartreeFock(newSystem)
{
    F = zeros<mat>(matDim, matDim);
    C = zeros<mat>(matDim, matDim);
    P = zeros<mat>(matDim, matDim);
    fockEnergy = ones<colvec>(matDim)*1.0E6;
}

//-------------------------------------------------------------------------------------------------------
// Solves the Hartree-Fock equations (iterated), stores the energy in double energy and the coefficients
// in vec C
void RHF::solve()
{
    double fockEnergyOld;
    double energyDiff = 1.0;

    // Calculate integrals
    C = zeros<mat>(matDim, nElectrons/2);
    calcIntegrals();


    // Iterate until the fock energy has converged
    while (energyDiff > toler){
        fockEnergyOld = fockEnergy(0);
        buildMatrix();
        solveSingle(F, C, P, fockEnergy);
        energyDiff = fabs(fockEnergyOld - fockEnergy(0));
    }

    // Calculate energy (not equal to Fock energy)
    energy = 0;

    for (int i = 0; i < matDim; i++){
        for (int j = 0; j < matDim; j++){
            energy += P(i,j)*h(i, j);
            for (int k = 0; k < matDim; k++){
                for (int l = 0; l < matDim; l++){
                    energy += 0.5*P(i,j)*P(l,k)*(Q[i][k][j][l] - 0.5*Q[i][k][l][j]);
                }
            }
        }
    }
    energy += system->getNucleiPotential();
}

//----------------------------------------------------------------------------------------------------------------
mat RHF::getCoeff(){
    return C;
}


//-----------------------------------------------------------------------------------------------------------------
// Builds F matrix (kinetic part, nuclear attraction part and electron-electron repulsion part, i.e.left hand side)

void RHF::buildMatrix()
{
    for (int i = 0; i < matDim; i++){
        for (int j = 0; j < matDim; j++){

            // One-electron integrals
            F(i,j) = h(i,j);

            // Add two-electron integrals
            for (int k = 0; k < matDim; k++){
                for (int l = 0; l < matDim; l++){
                    F(i,j) += 0.5*P(l,k)*(2*Q[i][k][j][l] - Q[i][k][l][j]);
                }
            }
        }
    }
}
