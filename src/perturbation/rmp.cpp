#include "rmp.h"

RMP::RMP(System* system, int perturbOrder, int frozenCore) :
    MollerPlesset(system, perturbOrder, frozenCore)
{
    m_solver = new RHF(m_system);
    m_nElectrons = m_system->getNumOfElectrons();
    m_fockEnergy = zeros(m_matDim);
    m_C = zeros(m_matDim, m_matDim);
    m_MOI.set_size(m_matDim, m_matDim);
    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            m_MOI(i,j) = zeros(m_matDim, m_matDim);
        }
    }
}

void RMP::solve()
{
    m_solver->solve();
    m_energyHF = m_solver->getEnergy();
    m_AOI = m_solver->getQmatrix();
    m_fockEnergy = m_solver->getFockEnergy()(0);
    m_C = m_solver->getCoeff()(0);

    // temp1, temp2, temp3 are temporary fields used to get from
    // Atomic Orbital Integrals to Molecular Orbital Integrals
    field<mat> temp1(m_matDim, m_matDim);
    field<mat> temp2(m_matDim, m_matDim);
    field<mat> temp3(m_matDim, m_matDim);
    // orbitalIntegral = <ij|g|ab>
    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            temp1(i,j) = zeros(m_matDim, m_matDim);
            temp2(i,j) = zeros(m_matDim, m_matDim);
            temp3(i,j) = zeros(m_matDim, m_matDim);
            m_MOI(i,j) = zeros(m_matDim, m_matDim);
        }
    }

    // Transform from Atomic Orbital Integrals to Molecular Orbital Integrals
    if (m_perturbOrder > 1){
//#ifdef RUN_MPI
//    clock_t begin = clock();
//#endif
        AOItoMOI(temp1, m_AOI, m_C, 0);
        AOItoMOI(temp2, temp1, m_C, 1);
        AOItoMOI(temp3, temp2, m_C, 2);
        AOItoMOI(m_MOI, temp3, m_C, 3);
//#ifdef RUN_MPI
//    clock_t end = clock();
//    int my_rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//    cout << "Proc " << my_rank << ": Time AOI to MOI: " << (double(end - begin))/CLOCKS_PER_SEC << endl;
//#endif
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

field<mat> RMP::getDensityMatrix()
{
    return m_solver->getDensityMatrix();
}


void RMP::calc2OrderPerturb()
{
    m_energy2Order = 0;
    // Sum up energy tems
    for (int i = m_frozenCore/2; i < m_nElectrons/2; i++){
        for (int j = m_frozenCore/2; j < m_nElectrons/2; j++){
            for (int a = m_nElectrons/2; a < m_matDim; a++){
                for (int b = m_nElectrons/2; b < m_matDim; b++){
                    m_energy2Order += m_MOI(i,j)(a,b)*(2*m_MOI(i,j)(a,b) - m_MOI(j,i)(a,b))
                                /(m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(a) - m_fockEnergy(b));
                }
            }
        }
    }
}


void RMP::calc3OrderPerturb()
{
    m_energy3Order = 0;
    // Contribution from particle ladder diagram
    for (int i = m_frozenCore/2; i < m_nElectrons/2; i++){
        for (int j = m_frozenCore/2; j < m_nElectrons/2; j++){
            for (int a = m_nElectrons/2; a < m_matDim; a++){
                for (int b = m_nElectrons/2; b < m_matDim; b++){
                    for (int c = m_nElectrons/2; c < m_matDim; c++){
                        for (int d = m_nElectrons/2; d < m_matDim; d++){
                            m_energy3Order += m_MOI(i,j)(a,b)*m_MOI(a,b)(c,d)*(2*m_MOI(c,d)(i,j) - m_MOI(c,d)(j,i))
                                          /((m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(a) - m_fockEnergy(b))
                                           *(m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(c) - m_fockEnergy(d)));
                        }
                    }
                }
            }
        }
    }

    // Contribution from hole ladder diagram
    for (int i = m_frozenCore/2; i < m_nElectrons/2; i++){
        for (int j = m_frozenCore/2; j < m_nElectrons/2; j++){
            for (int k = m_frozenCore/2; k < m_nElectrons/2; k++){
                for (int l = m_frozenCore/2; l < m_nElectrons/2; l++){
                    for (int a = m_nElectrons/2; a < m_matDim; a++){
                        for (int b = m_nElectrons/2; b < m_matDim; b++){
                            m_energy3Order += m_MOI(i,j)(a,b)*m_MOI(a,b)(k,l)*(2*m_MOI(k,l)(i,j) - m_MOI(k,l)(j,i))
                                           /((m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(a) - m_fockEnergy(b))
                                            *(m_fockEnergy(k) + m_fockEnergy(l) - m_fockEnergy(a) - m_fockEnergy(b)));
                        }
                    }
                }
            }
        }
    }


    // Contribution from loop diagram
    for (int i = m_frozenCore/2; i < m_nElectrons/2; i++){
        for (int j = m_frozenCore/2; j < m_nElectrons/2; j++){
            for (int k = m_frozenCore/2; k < m_nElectrons/2; k++){
                for (int a = m_nElectrons/2; a < m_matDim; a++){
                    for (int b = m_nElectrons/2; b < m_matDim; b++){
                        for (int c = m_nElectrons/2; c < m_matDim; c++){
                            m_energy3Order += -2*(m_MOI(i,j)(a,b)*m_MOI(k,b)(i,c)*(2*m_MOI(a,c)(k,j) - m_MOI(a,c)(j,k))
                                              +m_MOI(i,j)(a,b)*m_MOI(k,b)(c,i)*(2*m_MOI(a,c)(j,k) - m_MOI(a,c)(k,j))
                                              +m_MOI(i,j)(b,a)*m_MOI(k,b)(i,c)*(2*m_MOI(a,c)(j,k) - m_MOI(a,c)(k,j))
                                              +m_MOI(i,j)(b,a)*m_MOI(k,b)(c,i)*(2*m_MOI(a,c)(k,j) - 4*m_MOI(a,c)(j,k)))
                                            /((m_fockEnergy(i) + m_fockEnergy(j) - m_fockEnergy(a) - m_fockEnergy(b))
                                             *(m_fockEnergy(k) + m_fockEnergy(j) - m_fockEnergy(a) - m_fockEnergy(c)));
                        }
                    }
                }
            }
        }
    }
}
