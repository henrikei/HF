#include "hartreefock.h"

HartreeFock::HartreeFock(System *system):
    m_restrictedFactor(0)
{
    m_system = system;

    m_matDim = m_system->getTotalNumOfBasisFunc();
    m_h = zeros<mat>(m_matDim, m_matDim);
    m_S = zeros<mat>(m_matDim, m_matDim);
    m_Q.set_size(m_matDim, m_matDim);
    for (int p = 0; p < m_matDim; p++){
        for (int q = 0; q < m_matDim; q++){
            m_Q(p,q) = zeros(m_matDim, m_matDim);
        }
    }
    m_V = zeros(m_matDim, m_matDim);
    m_energy = 1.0E6;
    m_toler = 1.0E-8;
}


//----------------------------------------------------------------------------------------------------------------
double HartreeFock::getEnergy(){
    return m_energy;
}

field<mat> HartreeFock::getQmatrix()
{
    return m_Q;
}

//----------------------------------------------------------------------------------------------------------------
// Calculates all integrals needed to make the matrices
void HartreeFock::calcIntegrals()
{
    // The one-electron integrals (matrix h) and overlap inegrals (matrix S)
    rowvec oneElectronIntegrals;
    for (int i = 0; i < m_matDim; i++){
        for (int j = i; j < m_matDim; j++){
            oneElectronIntegrals = m_system->getOneElectronIntegrals(i,j);
            m_S(i,j) = oneElectronIntegrals(0);
            m_S(j,i) = m_S(i,j);
            m_h(i,j) = oneElectronIntegrals(1);
            m_h(j,i) = m_h(i,j);
        }
    }

    // The electron-electron integrals

    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < i+1; j++){
            for (int k = 0; k < i+1; k++){
                for (int l = 0; l < j+1; l++){
                    m_Q(i,j)(k,l) = m_system->getTwoElectronIntegral(i, j, k, l);
                    m_Q(k,j)(i,l) = m_Q(i,j)(k,l);
                    m_Q(i,l)(k,j) = m_Q(i,j)(k,l);
                    m_Q(k,l)(i,j) = m_Q(i,j)(k,l);

                    m_Q(j,i)(l,k) = m_Q(i,j)(k,l);
                    m_Q(j,k)(l,i) = m_Q(i,j)(k,l);
                    m_Q(l,i)(j,k) = m_Q(i,j)(k,l);
                    m_Q(l,k)(j,i) = m_Q(i,j)(k,l);
                }
            }
        }
    }
//        for (int i = 0; i < matDim; i++){
//            for (int j = 0; j < matDim; j++){
//                for (int k = 0; k < matDim; k++){
//                    for (int l = 0; l < matDim; l++){
//                        if ((Q[i][j][k][l] != Q[k][j][i][l])){
//                            cout << "Q not symmetric! 1" << endl;
//                        } else if (Q[i][j][k][l] != Q[i][l][k][j]){
//                            cout << "Q not symmetric! 2" << endl;
//                        } else if (Q[i][j][k][l] != Q[k][l][i][j]){
//                            cout << "Q not symmetric! 3" << endl;
//                        } else if(Q[i][j][k][l] != Q[j][i][l][k]){
//                            cout << "Q not symmetric! 4" << endl;
//                        } else if (Q[i][j][k][l] != Q[l][i][j][k]){
//                            cout << "Q not symmetric! 5" << endl;
//                        } else if (Q[i][j][k][l] != Q[j][k][l][i]){
//                            cout << "Q not symmetric! 6" << endl;
//                        } else if (Q[i][j][k][l] != Q[l][k][j][i]){
//                            cout << "Q not symmetric! 7" << endl;
//                        }
//                    }
//                }
//            }
    //        }
}

// Diagonalize overlap matrix S and calculate matrix m_V such that m_V.t()*S*m_V = I.
// Then set F2 = m_V.t()*F*m_V and C = m_V*C2.
void HartreeFock::diagOverlap()
{
    vec eigVal;

    eig_sym(eigVal, m_V, m_S);

    for (int i = 0; i < m_matDim; i++){
        m_V.col(i) = m_V.col(i)/sqrt(eigVal(i));
    }
}



//----------------------------------------------------------------------------------------------------------------
// Solves the Hartree-Fock equations (single iteration) and stores the Fock energy in double fockEnergy
// and coefficients in vec C;
void HartreeFock::solveSingle(const mat &Fock, mat &Coeffs, mat &P, colvec &fockEnergy, int nElectrons)
{
    vec eigVal;
    mat eigVec;
    mat F2 = zeros<mat>(m_matDim, m_matDim);

    F2 = m_V.t()*Fock*m_V;

    // Diagonalize matrix h2
    eig_sym(eigVal, eigVec, F2);
    Coeffs = m_V*eigVec;

    // Normalize the orbitals (phi = sum_p(C_p chi_p)). For vector C this means:
    double norm;
    for (int i = 0; i < m_matDim; i++){
        norm = dot(Coeffs.col(i), m_S*Coeffs.col(i));
        Coeffs.col(i) = Coeffs.col(i)/sqrt(norm);
    }

    // Compute density matrix. m_densityFactor accounts for the fact that the density matrix is
    // defined differently for RHF and UHF
    mat Ptemp = m_restrictedFactor*Coeffs.cols(0, nElectrons/m_restrictedFactor-1)*Coeffs.cols(0, nElectrons/m_restrictedFactor-1).t();
    P = 0.5*P + 0.5*Ptemp;  // Interpolate between new and old density matrix. Sometimes needed in order to achieve correct convergence.

    fockEnergy = eigVal;
}
