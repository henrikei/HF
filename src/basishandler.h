#ifndef BASISHANDLER_H
#define BASISHANDLER_H

#include <vector>
#include <armadillo>
#include "basisfunctions/basisfunctions.h"



class BasisHandler
{
public:
    BasisHandler();
    void addBasisFunctions(BasisFunctions* basis);

    int getTotalNumOfBasisFunc();
    rowvec3 getPosition(int p);
    rowvec getExponents(int p);
    rowvec getCoeffs(int p);
    irowvec getPowers(int p);
    int getAngMom(int p);
    int getAngMomMax();
private:
    vector<BasisFunctions*> allBasisFunctions;
    vector<int> map;
};

#endif // BASISHANDLER_H
