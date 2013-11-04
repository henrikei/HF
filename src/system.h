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
    rowvec2 getOneElectronIntegrals(int p, int q);
    double getTwoElectronIntegral(int p, int r, int q, int s);
    int getTotalNumOfBasisFunc();
private:
    int numberOfNuclei;
    mat nucleiPositions;
    BasisHandler* basisHandler;
    Integrator* integrator;
};

#endif // SYSTEM_H
