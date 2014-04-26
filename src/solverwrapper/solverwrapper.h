#ifndef SOLVERWRAPPER_H
#define SOLVERWRAPPER_H

#include "system/system.h"


template <typename T> class SolverWrapper
{
public:
    SolverWrapper(System *system);
    SolverWrapper(System *system, int perturbOrder);

    void solve();
    double getEnergy();
private:
    T *m_solver;
};



template <typename T>
SolverWrapper<T>::SolverWrapper(System *system)
{
    m_solver = new T(system);
}

template <typename T>
SolverWrapper<T>::SolverWrapper(System *system, int perturbOrder)
{
    m_solver = new T(system, perturbOrder);
}

template <typename T>
void SolverWrapper<T>::solve()
{
    m_solver->solve();
}

template <typename T>
double SolverWrapper<T>::getEnergy()
{
    return m_solver->getEnergy();
}

#endif // SOLVERWRAPPER_H
