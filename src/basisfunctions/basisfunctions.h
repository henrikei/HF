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
    rowvec getExponents(int p);
    rowvec getCoeffs(int p);
    rowvec3 getPowers(int p);
private:
    vector<rowvec> exponents;
    vector<rowvec> coeffs;
    vector<rowvec3> powers;
    rowvec3 position;
    int nBasisFunctions;
};

#endif // BASISFUNCTIONS_H
