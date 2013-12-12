#ifndef BASISFUNCTIONS_H
#define BASISFUNCTIONS_H

#include <vector>
#include <armadillo>
#include <boost/regex.hpp>

using namespace std;
using namespace arma;
using namespace boost;


class BasisFunctions
{
public:
    BasisFunctions();
    BasisFunctions(string inFileName);
    void setPosition(rowvec newPosition);
    rowvec3 getPosition();
    rowvec getExponents(int p);
    rowvec getCoeffs(int p);
    irowvec getPowers(int p);
    int getAngMom();
    int getNumOfBasisFunc();
protected:
    vector<rowvec> exponents;
    vector<rowvec> coeffs;
    vector<irowvec> powers;
    rowvec3 position;
    int nBasisFunctions;
    int angMom;

    void normalizeCoeffs();
    int factorial(int n);
};

#endif // BASISFUNCTIONS_H
