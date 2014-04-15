#include "mollerplesset.h"

MollerPlesset::MollerPlesset(System *system, int perturbOrder, int frozenCore)
{
    m_system = system;
    m_AOI.set_size(1,1);
    m_AOI(0,0) = zeros(1,1);
    m_perturbOrder = perturbOrder;
    m_matDim = m_system->getTotalNumOfBasisFunc();
    m_energyHF = 0;
    m_energy2Order = 0;
    m_energy3Order = 0;
    if (frozenCore % 2 == 0){
        m_frozenCore = frozenCore;
    } else {
        cout << "Error: The number of frozen core electrons must be even." << endl;
        exit(EXIT_FAILURE);
    }
}

double MollerPlesset::getEnergyHF()
{
    return m_energyHF;
}

double MollerPlesset::getEnergy2order()
{
    return m_energy2Order;
}

double MollerPlesset::getEnergy3order()
{
    return m_energy3Order;
}

// Transforms Atomic Orbital Integrals to Molecular Orbital Integrals one index at a time
void MollerPlesset::AOItoMOI(field<mat>& MOI, field<mat> AOI, mat C, int index)
{
    int a, b, c, d, e;

#ifdef RUN_MPI
    int numprocs, my_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    field<mat> MOI_local(m_matDim, m_matDim);
    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            MOI_local(i,j) = zeros(m_matDim, m_matDim);
        }
    }

    int N = m_matDim*m_matDim*m_matDim*m_matDim;
    int numcalcs = N/numprocs;
    int calcstart = my_rank*numcalcs;
    if (my_rank == numcalcs-1){
        numcalcs += N - numcalcs*numprocs;
    }
    int calcend = calcstart + numcalcs;

    int counter = 0;
    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            for (int k = 0; k < m_matDim; k++){
                for (int l = 0; l < m_matDim; l++){
                    if (calcstart-1 < counter && counter < calcend){
                        for (int m = 0; m < m_matDim; m++){
                            if (index == 0)       { e = i; a = m; b = j; c = k; d = l;}
                            else if(index == 1)   { e = j; a = i; b = m; c = k; d = l;}
                            else if(index == 2)   { e = k; a = i; b = j; c = m; d = l;}
                            else if(index == 3)   { e = l; a = i; b = j; c = k; d = m;}
                            MOI_local(i,j)(k,l) += C(m,e)*AOI(a,b)(c,d);
                        }
                    }
                    counter += 1;
                }
            }
        }
    }

    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            MPI_Allreduce(MOI_local(i,j).memptr(), MOI(i,j).memptr(), m_matDim*m_matDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }
    }
#endif

#ifndef RUN_MPI
    for (int i = 0; i < m_matDim; i++){
        for (int j = 0; j < m_matDim; j++){
            for (int k = 0; k < m_matDim; k++){
                for (int l = 0; l < m_matDim; l++){
                    for (int m = 0; m < m_matDim; m++){
                        if (index == 0)       { e = i; a = m; b = j; c = k; d = l;}
                        else if(index == 1)   { e = j; a = i; b = m; c = k; d = l;}
                        else if(index == 2)   { e = k; a = i; b = j; c = m; d = l;}
                        else if(index == 3)   { e = l; a = i; b = j; c = k; d = m;}
                        MOI(i,j)(k,l) += C(m,e)*AOI(a,b)(c,d);
                    }
                }
            }
        }
    }
#endif
}
