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

void BasisFunctions::normalizeCoeffs(){
    for (int i = 0; i < (int)coeffs.size(); i++){
        int l = sum(powers.at(i));
        for (int j = 0; j < (int)coeffs.at(i).n_elem; j++){
            if (l == 0){
                coeffs.at(i)(j) *= pow(2*exponents.at(i)(j)/M_PI, 0.75);
            } else if (l == 1){
                coeffs.at(i)(j) *= pow(2*exponents.at(i)(j)/M_PI, 0.75)*2*sqrt(exponents.at(i)(j));
            }
        }
    }
}
