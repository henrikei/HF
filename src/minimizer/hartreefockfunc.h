#ifndef HARTREEFOCKFUNC_H
#define HARTREEFOCKFUNC_H

#include <armadillo>
#include "func.h"
#include "hartreefock.h"
#include "system.h"

using namespace std;
using namespace arma;


class HartreeFockFunc : public Func
{
public:
    HartreeFockFunc(HartreeFock *solver, System *system);
    virtual rowvec getxInitial();
    virtual double getValue (rowvec x);
private:
    HartreeFock *m_solver;
    System *m_system;
};

#endif // HARTREEFOCKFUNC_H
