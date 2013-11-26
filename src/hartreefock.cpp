#include "hartreefock.h"

HartreeFock::HartreeFock(System *newSystem)
{
    system = newSystem;

    matDim = system->getTotalNumOfBasisFunc();
    nElectrons = system->getNumOfElectrons();
    h = zeros<mat>(matDim, matDim);
    F = zeros<mat>(matDim, matDim);
    S = zeros<mat>(matDim, matDim);
    C = zeros<mat>(matDim, nElectrons/2);
    P = zeros<mat>(matDim, matDim);
    fockEnergy = 1.0E6;
    energy = 1.0E6;
    toler = 1.0E-10;
}


//-------------------------------------------------------------------------------------------------------
// Solves the Hartree-Fock equations (iterated), stores the energy in double energy and the coefficients
// in vec C
void HartreeFock::solve()
{
    double fockEnergyOld;
    double energyDiff = 1.0;

    // Calculate integrals
    C = zeros<mat>(matDim, nElectrons/2);
    calcIntegrals();


    // Iterate until the fock energy has converged
    while (energyDiff > toler){
        fockEnergyOld = fockEnergy;
        buildMatrix();
        solveSingle();
        energyDiff = fabs(fockEnergyOld - fockEnergy);
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
double HartreeFock::getEnergy(){
    return energy;
}


//----------------------------------------------------------------------------------------------------------------
mat HartreeFock::getCoeff(){
    return C;
}


//-----------------------------------------------------------------------------------------------------------------
// Builds F matrix (kinetic part, nuclear attraction part and electron-electron repulsion part, i.e.left hand side)

void HartreeFock::buildMatrix()
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


//----------------------------------------------------------------------------------------------------------------
// Calculates all integrals needed to make the matrices
void HartreeFock::calcIntegrals()
{
    // The one-electron integrals (matrix h) and overlap inegrals (matrix S)
    rowvec oneElectronIntegrals;
    for (int i = 0; i < matDim; i++){
        for (int j = i; j < matDim; j++){
            oneElectronIntegrals = system->getOneElectronIntegrals(i,j);
            S(i,j) = oneElectronIntegrals(0);
            S(j,i) = S(i,j);
            h(i,j) = oneElectronIntegrals(1);
            h(j,i) = h(i,j);
        }
    }

    // The electron-electron integrals
    Q = new double***[matDim];
    for (int i = 0; i < matDim; i++){
        Q[i] = new double**[matDim];
        for (int j = 0; j < matDim; j++){
            Q[i][j] = new double*[matDim];
            for (int k = 0; k < matDim; k++){
                Q[i][j][k] = new double[matDim];
            }
        }
    }

    for (int i = 0; i < matDim; i++){
        for (int j = 0; j < i+1; j++){
            for (int k = 0; k < i+1; k++){
                for (int l = 0; l < j+1; l++){
                    Q[i][j][k][l] = system->getTwoElectronIntegral(i, j, k, l);
                    Q[k][j][i][l] = Q[i][j][k][l];
                    Q[i][l][k][j] = Q[i][j][k][l];
                    Q[k][l][i][j] = Q[i][j][k][l];

                    Q[j][i][l][k] = Q[i][j][k][l];
                    Q[j][k][l][i] = Q[i][j][k][l];
                    Q[l][i][j][k] = Q[i][j][k][l];
                    Q[l][k][j][i] = Q[i][j][k][l];
                }
            }
        }
    }
}



//----------------------------------------------------------------------------------------------------------------
// Solves the Hartree-Fock equations (single iteration) and stores the Fock energy in double fockEnergy
// and coefficients in vec C;
void HartreeFock::solveSingle()
{
    vec eigVal;
    mat eigVec;
    mat V = zeros<mat>(matDim, matDim);
    mat F2 = zeros<mat>(matDim, matDim);

    // Diagonalize overlap matrix S and calculate matrix V such that h2 = V.t()*h*V and C = V*C2

    eig_sym(eigVal, eigVec, S);

    for (int i = 0; i < matDim; i++){
        V.col(i) = eigVec.col(i)/sqrt(eigVal(i));
    }


    F2 = V.t()*F*V;

    // Diagonalize matrix h2

    eig_sym(eigVal, eigVec, F2);
    C = V*eigVec.cols(0, nElectrons/2-1);

    // Normalize vector C

    double norm;
    for (int i = 0; i < nElectrons/2; i++){
        norm = dot(C.col(i), S*C.col(i));
        C.col(i) = C.col(i)/sqrt(norm);
    }

    // Compute density matrix
    P = 2*C*C.t();

    fockEnergy = eigVal(0);
}
