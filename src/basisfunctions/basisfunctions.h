#ifndef BASISFUNCTIONS_H
#define BASISFUNCTIONS_H

#include <armadillo>
#include <boost/regex.hpp>
#include "primitive.h"
#include "contracted.h"

using namespace std;
using namespace arma;
using namespace boost;


class BasisFunctions
{
public:
    BasisFunctions();
    virtual ~BasisFunctions();
    void addContracteds(string inFileName, int posNum);
    void setPosPointer(mat* nucleiPositions);
    double getNumOfContracteds();
    Contracted* getContracted (int p);
    int getAngMomMax();

private:
    vector<Contracted*> m_contracteds;

    void addSomeContracteds(vector<double> exp, vector<double> coeff, field<imat> pows, field<rowvec> prefactor, int posNum);
    void normalizeCoeff(double exp, double& coeff, irowvec pows);
    int factorial (int n);
};

#endif // BASISFUNCTIONS_H
