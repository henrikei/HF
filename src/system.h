#ifndef SYSTEM_H
#define SYSTEM_H

#include <armadillo>
#include "integrator.h"
#include "basishandler.h"

using namespace std;
using namespace arma;



class System
{
public:
    System(BasisHandler* newBasisHandler, mat newNucleiPositions);
    double geth(int p, int q);
private:
    int numberOfNuclei;
    mat nucleiPositions;
    BasisHandler* basisHandler;
    Integrator* integrator;
};

#endif // SYSTEM_H
