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
    virtual rowvec getx();
    virtual double getValue (rowvec x);
    mat getNucleiPositions();
    void freezeAtoms(vector<int> frozenAtoms);

private:
    MollerPlesset *m_solver;
    System *m_system;

    mat m_trans;
    mat m_rotx;
    mat m_roty;
    mat m_rotz;

    mat m_nucleiPositions;
    mat m_nucleiPositionsTransformed;
    rowvec m_x;
    vector<pair<int,int>> m_map;

    void calcTransfMatrices();
    void transfPositions();
    void transfPositionsInverse();
    void matToVec();
    void vecToMat();
};

#endif // HARTREEFOCKFUNC_H
