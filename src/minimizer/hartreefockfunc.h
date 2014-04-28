#ifndef HARTREEFOCKFUNC_H
#define HARTREEFOCKFUNC_H

#include <armadillo>
#include "func.h"
#include "perturbation/mollerplesset.h"
#include "system/system.h"

using namespace std;
using namespace arma;


class HartreeFockFunc : public Func
{
public:
    HartreeFockFunc(MollerPlesset *solver, System *system);
    virtual rowvec getxInitial();
    virtual double getValue (rowvec x);
    mat getNucleiPositions();

private:
    MollerPlesset *m_solver;
    System *m_system;

    mat m_trans;
    mat m_rotx;
    mat m_roty;
    mat m_rotz;

    mat m_nucleiPositions;
    rowvec m_x;

    void calcTransfMatrices();
    void transfPositions();
    void transfPositionsInverse();
    void matToVec();
    void vecToMat();
};

#endif // HARTREEFOCKFUNC_H
