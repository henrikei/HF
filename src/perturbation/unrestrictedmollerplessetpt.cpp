#include "unrestrictedmollerplessetpt.h"

UnrestrictedMollerPlessetPT::UnrestrictedMollerPlessetPT(System *system, int perturbOrder, int frozenCore) :
    MollerPlessetPT(system, perturbOrder, frozenCore)
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

void UnrestrictedMollerPlessetPT::solve()
{
    m_solver->solve();
    m_energyHF = m_solver->getEnergy();
    m_AOI = m_solver->getQmatrix();
    m_fockEnergyUp = m_solver->getFockEnergy()(0);
    m_fockEnergyDown = m_solver->getFockEnergy()(1);
    m_Cup = m_solver->getCoeff()(0);
    m_Cdown = m_solver->getCoeff()(1);

    // Convert from Atomic Orbital Integrals (AOI) to Molecular Orbital Integrals (MOI)
    field<mat> tempUU1(m_matDim, m_matDim);
    field<mat> tempDD1(m_matDim, m_matDim);
    field<mat> tempUD1(m_matDim, m_matDim);
    field<mat> tempUU2(m_matDim, m_matDim);
    field<mat> tempDD2(m_matDim, m_matDim);
    field<mat> tempUD2(m_matDim, m_matDim);
    field<mat> tempUU3(m_matDim, m_matDim);
    field<mat> tempDD3(m_matDim, m_matDim);
    field<mat> tempUD3(m_matDim, m_matDim);

    for(int i = 0; i < m_matDim; i++){
        for(int j = 0; j < m_matDim; j++){
            tempUU1(i,j) = zeros(m_matDim, m_matDim);
            tempDD1(i,j) = zeros(m_matDim, m_matDim);
            tempUD1(i,j) = zeros(m_matDim, m_matDim);
            tempUU2(i,j) = zeros(m_matDim, m_matDim);
            tempDD2(i,j) = zeros(m_matDim, m_matDim);
            tempUD2(i,j) = zeros(m_matDim, m_matDim);
            tempUU3(i,j) = zeros(m_matDim, m_matDim);
            tempDD3(i,j) = zeros(m_matDim, m_matDim);
            tempUD3(i,j) = zeros(m_matDim, m_matDim);
        }
    }

    // Up-Up AOI to Up-Up MOI
    AOItoMOI(tempUU1, m_AOI, m_Cup, 0);
    AOItoMOI(tempUU2, tempUU1, m_Cup, 1);
    AOItoMOI(tempUU3, tempUU2, m_Cup, 2);
    AOItoMOI(m_MOI_UU, tempUU3, m_Cup, 3);

    // Down-Down AOI to Down-Down MOI
    AOItoMOI(tempDD1, m_AOI, m_Cdown, 0);
    AOItoMOI(tempDD2, tempDD1, m_Cdown, 1);
    AOItoMOI(tempDD3, tempDD2, m_Cdown, 2);
    AOItoMOI(m_MOI_DD, tempDD3, m_Cdown, 3);

    // Up-Down AOI to Up-Down MOI
    AOItoMOI(tempUD1, m_AOI, m_Cup, 0);
    AOItoMOI(tempUD2, tempUD1, m_Cdown, 1);
    AOItoMOI(tempUD3, tempUD2, m_Cup, 2);
    AOItoMOI(m_MOI_UD, tempUD3, m_Cdown, 3);

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

void UnrestrictedMollerPlessetPT::calc2OrderPerturb()
{
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

void UnrestrictedMollerPlessetPT::calc3OrderPerturb()
{
    // Contributions from loop diagram
    //-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
