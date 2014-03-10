#include "hartreefockfunc.h"

HartreeFockFunc::HartreeFockFunc(HartreeFock *solver, System *system)
{
    m_solver = solver;
    m_system = system;
}

rowvec HartreeFockFunc::getxInitial()
{
    mat nucleiPositions = m_system->getNucleiPositions();
}

double HartreeFockFunc::getValue(rowvec x)
{

}
