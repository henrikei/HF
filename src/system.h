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
    System(BasisHandler* newBasisHandler);
private:
    BasisHandler* basisHandler;
    mat h;
    mat F;
    mat S;
    vec C;
    field<mat> Q;

    double geth(int p, int q);
};

#endif // SYSTEM_H
