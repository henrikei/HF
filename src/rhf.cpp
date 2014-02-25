#include "rhf.h"

RHF::RHF(System* newSystem):HartreeFock(newSystem)
{
    F = zeros<mat>(matDim, matDim);
    C = zeros<mat>(matDim, matDim);
    P = zeros<mat>(matDim, matDim);
    fockEnergy = ones<colvec>(matDim)*1.0E6;
    perturbOrder = 1;
}

RHF::RHF(System* newSystem, int newPerturbOrder):HartreeFock(newSystem)
{
    F = zeros<mat>(matDim, matDim);
    C = zeros<mat>(matDim, matDim);
    P = zeros<mat>(matDim, matDim);
    fockEnergy = ones<colvec>(matDim)*1.0E6;
    perturbOrder = newPerturbOrder;
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

    // Perturbative terms
    if (perturbOrder == 2){
        energy += perturbation2order();
    } else if (perturbOrder == 1) {
    } else {
        cout << "Error. Only first and second order perturbation has been implemented." << endl;
        exit(EXIT_FAILURE);
    }
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


//--------------------------------------------------------------------------------------------------------------------
// Second order perturbation
double RHF::perturbation2order(){
    int nHStates = nElectrons/2;
    int nPStates = matDim - nElectrons/2;
    field<mat> orbitalIntegrals(nHStates, nHStates);               // orbitalIntegral = <ij|g|ab>
    for (int i = 0; i < nHStates; i++){
        for (int j = 0; j < nHStates; j++){
            orbitalIntegrals(i,j) = zeros(nPStates, nPStates);
            for (int a = 0; a < nPStates; a++){
                for (int b = 0; b < nPStates; b++){
                    for (int p = 0; p < matDim; p++){
                        for (int q = 0; q < matDim; q++){
                            for (int r = 0; r < matDim; r++){
                                for (int s = 0; s < matDim; s++){
                                    orbitalIntegrals(i,j)(a,b) += C(p,i)*C(q,j)*C(r,a+nHStates)*C(s,b+nHStates)*Q[p][q][r][s];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < nHStates; i++){
        for (int j = 0; j < nHStates; j++){
            for (int a = 0; a < nPStates; a++){
                for (int b = 0; b < nPStates; b++){
                    energyMP2 += orbitalIntegrals(i,j)(a,b)*(2*orbitalIntegrals(i,j)(a,b) - orbitalIntegrals(j,i)(a,b))
                                /(fockEnergy(i) + fockEnergy(j) - fockEnergy(a+nHStates) - fockEnergy(b+nHStates));
                }
            }
        }
    }

    return energyMP2;

}
