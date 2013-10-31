#ifndef SYSTEM_H
#define SYSTEM_H

#include <armadillo>
#include "integrator.h"

using namespace std;
using namespace arma;



class System
{
public:
    System(const mat &newNucleiPos, const vector<rowvec> &newExponents, const vector<rowvec> &newCoeffs);

private:
    mat nucleiPos;
    vector<rowvec> exponents;
    vector<rowvec> coeffs;
    int nNuclei;
    int nBasisFunc;
    int matDim;

    Integrator integrator;

    double getOneParticleInt(int p, int q);
    double getTwoParticleInt(int p, int q, int r, int s);
};

#endif // SYSTEM_H
