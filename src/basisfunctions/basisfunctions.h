#ifndef BASISFUNCTIONS_H
#define BASISFUNCTIONS_H

#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;



class BasisFunctions
{
public:
    BasisFunctions();
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
};

#endif // BASISFUNCTIONS_H
