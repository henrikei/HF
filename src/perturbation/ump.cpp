#include "ump.h"

UMP::UMP(System *system, int perturbOrder, int frozenCore) :
    MollerPlesset(system, perturbOrder, frozenCore)
{
    m_solver = new UHF(m_system);
    m_nElectronsUp = m_solver->getNumOfElectronsUp();
    m_nElectronsDown = m_solver->getNumOfElecrtonsDown();
    m_fockEnergyUp = zeros(m_matDim);
    m_fockEnergyDown = zeros(m_matDim);
    m_Cup = zeros(m_matDim, m_matDim);
    m_Cdown = zeros(m_matDim, m_matDim);
    m_MOI_UU.set_size(m_matDim, m_matDim);
    m_MOI_DD.set_size(m_matDim, m_matDim);
    m_MOI_UD.set_size(m_matDim, m_matDim);
    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            m_MOI_UU(i,j) = zeros(m_matDim, m_matDim);
            m_MOI_DD(i,j) = zeros(m_matDim, m_matDim);
            m_MOI_UD(i,j) = zeros(m_matDim, m_matDim);
        }
    }
}

void UMP::solve()
{
    m_solver->solve();
    m_energyHF = m_solver->getEnergy();
    m_AOI = m_solver->getQmatrix();
    m_fockEnergyUp = m_solver->getFockEnergy()(0);
    m_fockEnergyDown = m_solver->getFockEnergy()(1);
    m_Cup = m_solver->getCoeff()(0);
    m_Cdown = m_solver->getCoeff()(1);

    // Convert from Atomic Orbital Integrals (AOI) to Molecular Orbital Integrals (MOI)
    if (m_perturbOrder > 1){
        field<mat> temp1(m_matDim, m_matDim);
        field<mat> temp2(m_matDim, m_matDim);
        cout << "Done initializing temps" << endl;

        for(int i = 0; i < m_matDim; i++){
            for(int j = 0; j < m_matDim; j++){
                temp1(i,j) = zeros(m_matDim, m_matDim);
                temp2(i,j) = zeros(m_matDim, m_matDim);
                m_MOI_UU(i,j) = zeros(m_matDim, m_matDim);
                m_MOI_DD(i,j) = zeros(m_matDim, m_matDim);
                m_MOI_UD(i,j) = zeros(m_matDim, m_matDim);
            }
        }
        cout << "Done initializing temps" << endl;

        // Up-Up AOI to Up-Up MOI
        AOItoMOI(temp1, m_AOI, m_Cup, 0);
        AOItoMOI(temp2, temp1, m_Cup, 1);
        fillZero(temp1);
        AOItoMOI(temp1, temp2, m_Cup, 2);
        AOItoMOI(m_MOI_UU, temp1, m_Cup, 3);
        fillZero(temp1);
        fillZero(temp2);
        cout << "1" << endl;

        // Down-Down AOI to Down-Down MOI
        AOItoMOI(temp1, m_AOI, m_Cdown, 0);
        AOItoMOI(temp2, temp1, m_Cdown, 1);
        fillZero(temp1);
        AOItoMOI(temp1, temp2, m_Cdown, 2);
        AOItoMOI(m_MOI_DD, temp1, m_Cdown, 3);
        fillZero(temp1);
        fillZero(temp2);
        cout << "2" << endl;

        // Up-Down AOI to Up-Down MOI
        AOItoMOI(temp1, m_AOI, m_Cup, 0);
        AOItoMOI(temp2, temp1, m_Cdown, 1);
        fillZero(temp1);
        AOItoMOI(temp1, temp2, m_Cup, 2);
        AOItoMOI(m_MOI_UD, temp1, m_Cdown, 3);
        cout << "3" << endl;
    }

    if (m_perturbOrder == 1){
    } else if (m_perturbOrder == 2){
        calc2OrderPerturb();
    } else if (m_perturbOrder == 3){
        calc2OrderPerturb();
        calc3OrderPerturb();
    } else {
        cout << "Error: " << m_perturbOrder << " order perturbation theory not implemented." << endl;
        exit(EXIT_FAILURE);
    }
}

field<mat> UMP::getDensityMatrix()
{
    return m_solver->getDensityMatrix();
}

void UMP::calc2OrderPerturb()
{
    m_energy2Order = 0;
    // Sum up energy terms
    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsUp; j++){
            for (int a = m_nElectronsUp; a < m_matDim; a++){
                for (int b = m_nElectronsUp; b < m_matDim; b++){
                    m_energy2Order += 0.25*(m_MOI_UU(i,j)(a,b) - m_MOI_UU(i,j)(b,a))*(m_MOI_UU(a,b)(i,j) - m_MOI_UU(b,a)(i,j))/
                            (m_fockEnergyUp(i) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyUp(b));
                }
            }
        }
    }
    for (int i = m_frozenCore/2; i < m_nElectronsDown; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int a = m_nElectronsDown; a < m_matDim; a++){
                for (int b = m_nElectronsDown; b < m_matDim; b++){
                    m_energy2Order += 0.25*(m_MOI_DD(i,j)(a,b) - m_MOI_DD(i,j)(b,a))*(m_MOI_DD(a,b)(i,j) - m_MOI_DD(b,a)(i,j))/
                            (m_fockEnergyDown(i) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyDown(b));
                }
            }
        }
    }
    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int a = m_nElectronsUp; a < m_matDim; a++){
                for (int b = m_nElectronsDown; b < m_matDim; b++){
                    m_energy2Order += m_MOI_UD(i,j)(a,b)*m_MOI_UD(i,j)(a,b)/
                            (m_fockEnergyUp(i) + m_fockEnergyDown(j) - m_fockEnergyUp(a) - m_fockEnergyDown(b));
                }
            }
        }
    }
}

void UMP::calc3OrderPerturb()
{
    // Contributions from loop diagram
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    m_energy3Order = 0;
    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsUp; j++){
            for (int k = m_frozenCore/2; k < m_nElectronsUp; k++){
                for (int a = m_nElectronsUp; a < m_matDim; a++){
                    for (int b = m_nElectronsUp; b < m_matDim; b++){
                        for (int c = m_nElectronsUp; c < m_matDim; c++){
                            m_energy3Order += -(m_MOI_UU(i,j)(a,b) - m_MOI_UU(i,j)(b,a))*(m_MOI_UU(k,b)(i,c) - m_MOI_UU(k,b)(c,i))*(m_MOI_UU(a,c)(k,j) - m_MOI_UU(a,c)(j,k))
                                          /((m_fockEnergyUp(i) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyUp(b))
                                            *(m_fockEnergyUp(k) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyUp(c)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int k = m_frozenCore/2; k < m_nElectronsUp; k++){
                for (int a = m_nElectronsDown; a < m_matDim; a++){
                    for (int b = m_nElectronsUp; b < m_matDim; b++){
                        for (int c = m_nElectronsUp; c < m_matDim; c++){
                            m_energy3Order += -m_MOI_UD(i,j)(b,a)*(m_MOI_UU(k,b)(i,c) - m_MOI_UU(k,b)(c,i))*m_MOI_UD(c,a)(k,j)
                                          /((m_fockEnergyUp(i) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyUp(b))
                                            *(m_fockEnergyUp(k) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyUp(c)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsDown; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsUp; j++){
            for (int k = m_frozenCore/2; k < m_nElectronsUp; k++){
                for (int a = m_nElectronsUp; a < m_matDim; a++){
                    for (int b = m_nElectronsDown; b < m_matDim; b++){
                        for (int c = m_nElectronsUp; c < m_matDim; c++){
                            m_energy3Order += -m_MOI_UD(j,i)(a,b)*m_MOI_UD(k,b)(c,i)*(m_MOI_UU(a,c)(k,j) - m_MOI_UU(a,c)(j,k))
                                          /((m_fockEnergyDown(i) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyDown(b))
                                            *(m_fockEnergyUp(k) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyUp(c)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsDown; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int k = m_frozenCore/2; k < m_nElectronsUp; k++){
                for (int a = m_nElectronsDown; a < m_matDim; a++){
                    for (int b = m_nElectronsDown; b < m_matDim; b++){
                        for (int c = m_nElectronsUp; c < m_matDim; c++){
                            m_energy3Order += -(m_MOI_DD(i,j)(a,b) - m_MOI_DD(i,j)(b,a))*m_MOI_UD(k,b)(c,i)*m_MOI_UD(c,a)(k,j)
                                          /((m_fockEnergyDown(i) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyDown(b))
                                            *(m_fockEnergyUp(k) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyUp(c)));
                        }
                    }
                }
            }
        }
    }

    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int k = m_frozenCore/2; k < m_nElectronsUp; k++){
                for (int a = m_nElectronsUp; a < m_matDim; a++){
                    for (int b = m_nElectronsDown; b < m_matDim; b++){
                        for (int c = m_nElectronsDown; c < m_matDim; c++){
                            m_energy3Order += -m_MOI_UD(i,j)(a,b)*m_MOI_UD(k,b)(i,c)*m_MOI_UD(a,c)(k,j)
                                          /((m_fockEnergyUp(i) + m_fockEnergyDown(j) - m_fockEnergyUp(a) - m_fockEnergyDown(b))
                                            *(m_fockEnergyUp(k) + m_fockEnergyDown(j) - m_fockEnergyUp(a) - m_fockEnergyDown(c)));
                        }
                    }
                }
            }
        }
    }



    for (int i = m_frozenCore/2; i < m_nElectronsDown; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int k = m_frozenCore/2; k < m_nElectronsDown; k++){
                for (int a = m_nElectronsDown; a < m_matDim; a++){
                    for (int b = m_nElectronsDown; b < m_matDim; b++){
                        for (int c = m_nElectronsDown; c < m_matDim; c++){
                            m_energy3Order += -(m_MOI_DD(i,j)(a,b) - m_MOI_DD(i,j)(b,a))*(m_MOI_DD(k,b)(i,c) - m_MOI_DD(k,b)(c,i))*(m_MOI_DD(a,c)(k,j) - m_MOI_DD(a,c)(j,k))
                                          /((m_fockEnergyDown(i) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyDown(b))
                                            *(m_fockEnergyDown(k) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyDown(c)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsDown; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsUp; j++){
            for (int k = m_frozenCore/2; k < m_nElectronsDown; k++){
                for (int a = m_nElectronsUp; a < m_matDim; a++){
                    for (int b = m_nElectronsDown; b < m_matDim; b++){
                        for (int c = m_nElectronsDown; c < m_matDim; c++){
                            m_energy3Order += -m_MOI_UD(j,i)(a,b)*(m_MOI_DD(k,b)(i,c) - m_MOI_DD(k,b)(c,i))*m_MOI_UD(a,c)(j,k)
                                          /((m_fockEnergyDown(i) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyDown(b))
                                            *(m_fockEnergyDown(k) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyDown(c)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int k = m_frozenCore/2; k < m_nElectronsDown; k++){
                for (int a = m_nElectronsDown; a < m_matDim; a++){
                    for (int b = m_nElectronsUp; b < m_matDim; b++){
                        for (int c = m_nElectronsDown; c < m_matDim; c++){
                            m_energy3Order += -m_MOI_UD(i,j)(b,a)*m_MOI_UD(b,k)(i,c)*(m_MOI_DD(a,c)(k,j) - m_MOI_DD(a,c)(j,k))
                                          /((m_fockEnergyUp(i) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyUp(b))
                                            *(m_fockEnergyDown(k) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyDown(c)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsUp; j++){
            for (int k = m_frozenCore/2; k < m_nElectronsDown; k++){
                for (int a = m_nElectronsUp; a < m_matDim; a++){
                    for (int b = m_nElectronsUp; b < m_matDim; b++){
                        for (int c = m_nElectronsDown; c < m_matDim; c++){
                            m_energy3Order += -(m_MOI_UU(i,j)(a,b) - m_MOI_UU(i,j)(b,a))*m_MOI_UD(b,k)(i,c)*m_MOI_UD(a,c)(j,k)
                                          /((m_fockEnergyUp(i) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyUp(b))
                                            *(m_fockEnergyDown(k) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyDown(c)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsDown; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsUp; j++){
            for (int k = m_frozenCore/2; k < m_nElectronsDown; k++){
                for (int a = m_nElectronsDown; a < m_matDim; a++){
                    for (int b = m_nElectronsUp; b < m_matDim; b++){
                        for (int c = m_nElectronsUp; c < m_matDim; c++){
                            m_energy3Order += -m_MOI_UD(j,i)(b,a)*m_MOI_UD(b,k)(c,i)*m_MOI_UD(c,a)(j,k)
                                          /((m_fockEnergyDown(i) + m_fockEnergyUp(j) - m_fockEnergyDown(a) - m_fockEnergyUp(b))
                                            *(m_fockEnergyDown(k) + m_fockEnergyUp(j) - m_fockEnergyDown(a) - m_fockEnergyUp(c)));
                        }
                    }
                }
            }
        }
    }



    // Contributions from particle ladder diagram
    // -----------------------------------------------------------------------------------------------------------------------------------------------------------------

    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsUp; j++){
            for (int a = m_nElectronsUp; a < m_matDim; a++){
                for (int b = m_nElectronsUp; b < m_matDim; b++){
                    for (int c = m_nElectronsUp; c < m_matDim; c++){
                        for (int d = m_nElectronsUp; d < m_matDim; d++){
                            m_energy3Order += (m_MOI_UU(i,j)(a,b) - m_MOI_UU(i,j)(b,a))*(m_MOI_UU(a,b)(c,d) - m_MOI_UU(a,b)(d,c))*(m_MOI_UU(c,d)(i,j) - m_MOI_UU(c,d)(j,i))
                                           /(8*(m_fockEnergyUp(i) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyUp(b))
                                            *(m_fockEnergyUp(i) + m_fockEnergyUp(j) - m_fockEnergyUp(c) - m_fockEnergyUp(d)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsDown; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int a = m_nElectronsDown; a < m_matDim; a++){
                for (int b = m_nElectronsDown; b < m_matDim; b++){
                    for (int c = m_nElectronsDown; c < m_matDim; c++){
                        for (int d = m_nElectronsDown; d < m_matDim; d++){
                            m_energy3Order += (m_MOI_DD(i,j)(a,b) - m_MOI_DD(i,j)(b,a))*(m_MOI_DD(a,b)(c,d) - m_MOI_DD(a,b)(d,c))*(m_MOI_DD(c,d)(i,j) - m_MOI_DD(c,d)(j,i))
                                           /(8*(m_fockEnergyDown(i) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyDown(b))
                                            *(m_fockEnergyDown(i) + m_fockEnergyDown(j) - m_fockEnergyDown(c) - m_fockEnergyDown(d)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int a = m_nElectronsUp; a < m_matDim; a++){
                for (int b = m_nElectronsDown; b < m_matDim; b++){
                    for (int c = m_nElectronsUp; c < m_matDim; c++){
                        for (int d = m_nElectronsDown; d < m_matDim; d++){
                            m_energy3Order += m_MOI_UD(i,j)(a,b)*m_MOI_UD(a,b)(c,d)*m_MOI_UD(c,d)(i,j)
                                           /((m_fockEnergyUp(i) + m_fockEnergyDown(j) - m_fockEnergyUp(a) - m_fockEnergyDown(b))
                                            *(m_fockEnergyUp(i) + m_fockEnergyDown(j) - m_fockEnergyUp(c) - m_fockEnergyDown(d)));
                        }
                    }
                }
            }
        }
    }



    // Contributions from hole ladder diagram
    // --------------------------------------------------------------------------------------------------------------------------------------------------------------------

    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsUp; j++){
            for (int a = m_nElectronsUp; a < m_matDim; a++){
                for (int b = m_nElectronsUp; b < m_matDim; b++){
                    for (int k = m_frozenCore/2; k < m_nElectronsUp; k++){
                        for (int l = m_frozenCore/2; l < m_nElectronsUp; l++){
                            m_energy3Order += (m_MOI_UU(i,j)(a,b) - m_MOI_UU(i,j)(b,a))*(m_MOI_UU(a,b)(k,l) - m_MOI_UU(a,b)(l,k))*(m_MOI_UU(k,l)(i,j) - m_MOI_UU(k,l)(j,i))
                                           /(8*(m_fockEnergyUp(i) + m_fockEnergyUp(j) - m_fockEnergyUp(a) - m_fockEnergyUp(b))
                                            *(m_fockEnergyUp(k) + m_fockEnergyUp(l) - m_fockEnergyUp(a) - m_fockEnergyUp(b)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsDown; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int a = m_nElectronsDown; a < m_matDim; a++){
                for (int b = m_nElectronsDown; b < m_matDim; b++){
                    for (int k = m_frozenCore/2; k < m_nElectronsDown; k++){
                        for (int l = m_frozenCore/2; l < m_nElectronsDown; l++){
                            m_energy3Order += (m_MOI_DD(i,j)(a,b) - m_MOI_DD(i,j)(b,a))*(m_MOI_DD(a,b)(k,l) - m_MOI_DD(a,b)(l,k))*(m_MOI_DD(k,l)(i,j) - m_MOI_DD(k,l)(j,i))
                                           /(8*(m_fockEnergyDown(i) + m_fockEnergyDown(j) - m_fockEnergyDown(a) - m_fockEnergyDown(b))
                                            *(m_fockEnergyDown(k) + m_fockEnergyDown(l) - m_fockEnergyDown(a) - m_fockEnergyDown(b)));
                        }
                    }
                }
            }
        }
    }


    for (int i = m_frozenCore/2; i < m_nElectronsUp; i++){
        for (int j = m_frozenCore/2; j < m_nElectronsDown; j++){
            for (int a = m_nElectronsUp; a < m_matDim; a++){
                for (int b = m_nElectronsDown; b < m_matDim; b++){
                    for (int k = m_frozenCore/2; k < m_nElectronsUp; k++){
                        for (int l = m_frozenCore/2; l < m_nElectronsDown; l++){
                            m_energy3Order += m_MOI_UD(i,j)(a,b)*m_MOI_UD(a,b)(k,l)*m_MOI_UD(k,l)(i,j)
                                           /((m_fockEnergyUp(i) + m_fockEnergyDown(j) - m_fockEnergyUp(a) - m_fockEnergyDown(b))
                                            *(m_fockEnergyUp(k) + m_fockEnergyDown(l) - m_fockEnergyUp(a) - m_fockEnergyDown(b)));
                        }
                    }
                }
            }
        }
    }
}
