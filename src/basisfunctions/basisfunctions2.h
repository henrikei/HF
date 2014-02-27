#ifndef BASISFUNCTIONS2_H
#define BASISFUNCTIONS2_H

#include <armadillo>
#include <boost/regex.hpp>
#include "primitive.h"
#include "contracted.h"

using namespace std;
using namespace arma;
using namespace boost;


class BasisFunctions2
{
public:
    BasisFunctions2();
    void addContracteds(string inFileName, rowvec3 position);
    double getNumOfContracteds();
    Contracted* getContracted (int p);
    int getAngMomMax();
private:
    vector<Contracted*> m_contracteds;

    void normalizeCoeff(double exp, double& coeff, irowvec pows);
    int factorial (int n);
};

#endif // BASISFUNCTIONS2_H
