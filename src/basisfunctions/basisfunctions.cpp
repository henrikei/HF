#include "basisfunctions.h"


BasisFunctions::BasisFunctions()
{
}


void BasisFunctions::setPosition(rowvec newPosition)
{
    position = newPosition;
}


rowvec BasisFunctions::getExponents(int p)
{
    return exponents.at(p);
}


rowvec BasisFunctions::getCoeffs(int p)
{
    return coeffs.at(p);
}

rowvec3 BasisFunctions::getPowers(int p)
{
    return powers.at(p);
}

int BasisFunctions::getNumOfBasisFunc()
{
    return exponents.size();
}
