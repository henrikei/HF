#include "system.h"



System::System(BasisHandler *newBasisHandler, mat newNucleiPositions)
{
    basisHandler = newBasisHandler;
    integrator = new Integrator;

    nucleiPositions = newNucleiPositions;

    matDim = basisHandler->getTotalNumOfBasisFunc();
}


double System::geth(int p, int q)
{
    rowvec expA = basisHandler->getExponents(p);
    rowvec expB = basisHandler->getExponents(q);

    rowvec coeffsA = basisHandler->getCoeffs(p);
    rowvec coeffsB = basisHandler->getCoeffs(q);

    irowvec powersA = basisHandler->getPowers(p);
    irowvec powersB = basisHandler->getPowers(q);

    rowvec3 positionA = basisHandler->getPosition(p);
    rowvec3 positionB = basisHandler->getPosition(q);

    int angMom = basisHandler->getAngMom(p);

    int nPrimitives = expA.n_elem;


    integrator->setPositionA(positionA);
    integrator->setPositionB(positionB);

    integrator->setMaxAngMom(angMom);

    int i = powersA(0);
    int j = powersB(0);
    int k = powersA(1);
    int l = powersB(1);
    int m = powersA(2);
    int n = powersB(2);

    double value = 0;

    for (int v = 0; v < nPrimitives ; v++){
        for (int w = 0; w < nPrimitives; w++){
            integrator->setAlpha(expA(v));
            integrator->setBeta(expB(w));
            integrator->setE();
            value += integrator->kinetic(i, j, k, l, m, n)*coeffsA(v)*coeffsB(w);
        }
    }

    return value;
}
