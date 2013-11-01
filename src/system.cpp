#include "system.h"



System::System(BasisHandler *newBasisHandler)
{
    basisHandler = newBasisHandler;

    int matDim = basisHandler->getTotalNumOfBasisFunc();

    h = zeros<mat>(matDim, matDim);
    F = zeros<mat>(matDim, matDim);
    S = zeros<mat>(matDim, matDim);
    C = zeros<vec>(matDim);

    Q.set_size(matDim, matDim);
    for(int i = 0; i < matDim; i++){
        for(int j = 0; j < matDim; j++){
            Q(i,j) = zeros<mat>(matDim, matDim);
        }
    }
}

double System::geth(int p, int q)
{
    rowvec expA = basisHandler->getExponents(p);
    rowvec ecpB = basisHandler->getExponents(q);

    rowvec coeffsA = basisHandler->getCoeffs(p);
    rowvec coeffsB = basisHandler->getCoeffs(q);

    irowvec powersA = basisHandler->getPowers(p);
    irowvec powersB = basisHandler->getPowers(q);

    int n = expA.n_elem;
}
