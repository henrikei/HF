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
            } else if (l == 2){
                int lx = powers.at(i)(0);
                int ly = powers.at(i)(1);
                int lz = powers.at(i)(2);
                coeffs.at(i)(j) *= pow(2*exponents.at(i)(j)/M_PI, 0.75)*sqrt(pow(8*exponents.at(i)(j),lx+ly+lz)
                                   *factorial(2*lx)*factorial(2*ly)*factorial(2*lz)/(factorial(2*lx)*factorial(2*ly)*factorial(2*lz)));
            }
        }
    }
}

int BasisFunctions::factorial(int n)
{
    double value = 1;
    double i = 1;

    while(i < n){
        i += 1;
        value *= i;
    }
    return value;
}
