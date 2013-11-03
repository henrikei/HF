#include "basisfunctions.h"


BasisFunctions::BasisFunctions()
{
}


void BasisFunctions::setPosition(rowvec newPosition)
{
    position = newPosition;
}

rowvec3 BasisFunctions::getPosition()
{
    return position;
}


rowvec BasisFunctions::getExponents(int p)
{
    return exponents.at(p);
}


rowvec BasisFunctions::getCoeffs(int p)
{
    return coeffs.at(p);
}

irowvec BasisFunctions::getPowers(int p)
{
    return powers.at(p);
}

int BasisFunctions::getAngMom(){
    return angMom;
}

int BasisFunctions::getNumOfBasisFunc()
{
    return exponents.size();
}
