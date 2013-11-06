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
    System(BasisHandler* newBasisHandler, mat newNucleiPositions, rowvec newCharges, int newNElectrons);
    rowvec2 getOneElectronIntegrals(int p, int q);
    double getTwoElectronIntegral(int p, int q, int r, int s);
    double getNucleiPotential();
    int getTotalNumOfBasisFunc();
    int getNumOfElectrons();
private:
    int numberOfNuclei;
    mat nucleiPositions;
    rowvec charges;
    int nElectrons;
    BasisHandler* basisHandler;
    Integrator* integrator;
};

#endif // SYSTEM_H
