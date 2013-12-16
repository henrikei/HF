#include "hartreefock.h"

HartreeFock::HartreeFock(System *newSystem)
{
    system = newSystem;

    matDim = system->getTotalNumOfBasisFunc();
    nElectrons = system->getNumOfElectrons();
    h = zeros<mat>(matDim, matDim);
    S = zeros<mat>(matDim, matDim);
    energy = 1.0E6;
    energyMP2 = 0.0;
    toler = 1.0E-10;
}

//----------------------------------------------------------------------------------------------------------------
double HartreeFock::getEnergy(){
    return energy;
}

double HartreeFock::getEnergyMP2()
{
    return energyMP2;
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
        for (int j = 0; j < matDim; j++){
            for (int k = 0; k < matDim; k++){
                for (int l = 0; l < matDim; l++){
                    Q[i][j][k][l] = system->getTwoElectronIntegral(i, j, k, l);
                }
            }
        }
    }

//    for (int i = 0; i < matDim; i++){
//        for (int j = 0; j < matDim; j++){
//            for (int k = 0; k < matDim; k++){
//                for (int l = 0; l < matDim; l++){
//                    if (fabs(Q[i][j][k][l] - Q[k][j][i][l]) > 1.0E-8){
//                        cout << "Q not symmetric! 1" << endl;
//                    } else if (fabs(Q[i][j][k][l] - Q[i][l][k][j]) > 1.0E-8){
//                        cout << "Q not symmetric! 2" << endl;
//                    } else if (fabs(Q[i][j][k][l] - Q[k][l][i][j]) > 1.0E-8){
//                        cout << "Q not symmetric! 2" << endl;
//                    } else if (fabs(Q[i][j][k][l] - Q[j][i][l][k]) > 1.0E-8){
//                        cout << "Q not symmetric! 2" << endl;
//                    }  else if (fabs(Q[i][j][k][l] - Q[l][i][j][k]) > 1.0E-8){
//                        cout << "Q not symmetric! 2" << endl;
//                    }  else if (fabs(Q[i][j][k][l] - Q[j][k][l][i]) > 1.0E-8){
//                        cout << "Q not symmetric! 2" << endl;
//                    }  else if (fabs(Q[i][j][k][l] - Q[l][k][j][i]) > 1.0E-8){
//                        cout << "Q not symmetric! 2" << endl;
//                    }
//                }
//            }
//        }
//    }

//    for (int i = 0; i < matDim; i++){
//        for (int j = 0; j < i+1; j++){
//            for (int k = 0; k < i+1; k++){
//                for (int l = 0; l < j+1; l++){
//                    Q[i][j][k][l] = system->getTwoElectronIntegral(i, j, k, l);
//                    Q[k][j][i][l] = Q[i][j][k][l];
//                    Q[i][l][k][j] = Q[i][j][k][l];
//                    Q[k][l][i][j] = Q[i][j][k][l];

//                    Q[j][i][l][k] = Q[i][j][k][l];
//                    Q[j][k][l][i] = Q[i][j][k][l];
//                    Q[l][i][j][k] = Q[i][j][k][l];
//                    Q[l][k][j][i] = Q[i][j][k][l];
//                }
//            }
//        }
//    }
}



//----------------------------------------------------------------------------------------------------------------
// Solves the Hartree-Fock equations (single iteration) and stores the Fock energy in double fockEnergy
// and coefficients in vec C;
void HartreeFock::solveSingle(const mat &Fock, mat &Coeffs, mat &P, colvec &fockEnergy)
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

    F2 = V.t()*Fock*V;

    // Diagonalize matrix h2
    eig_sym(eigVal, eigVec, F2);
    Coeffs = V*eigVec;

    // Normalize vector C
    double norm;
    for (int i = 0; i < matDim; i++){
        norm = dot(Coeffs.col(i), S*Coeffs.col(i));
        Coeffs.col(i) = Coeffs.col(i)/sqrt(norm);
    }

    // Compute density matrix
    mat Ptemp = 2*Coeffs.cols(0, nElectrons/2-1)*Coeffs.cols(0, nElectrons/2-1).t();
    P = 0.5*P + 0.5*Ptemp;  // Interpolate between new and old density matrix

    fockEnergy = eigVal;
}
