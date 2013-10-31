#include "system.h"



System::System(const mat &newNucleiPos, const vector<rowvec> &newExponents, const vector<rowvec> &newCoeffs)
{
    nucleiPos = newNucleiPos;
    exponents = newExponents;
    coeffs = newCoeffs;

    nNuclei = nucleiPos.n_rows;
    nBasisFunc = exponents.size();
    matDim = nNuclei*nBasisFunc;
}

double System::getOneParticleInt(int p, int q)
{
    int p1 = p/nNuclei;
    int p2 = p%nNuclei;
    int q1 = q/nNuclei;
    int q2 = q/nNuclei;


}

double System::getTwoParticleInt(int p, int q, int r, int s)
{
}


