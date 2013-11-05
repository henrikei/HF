#include "system.h"



System::System(BasisHandler *newBasisHandler, mat newNucleiPositions, rowvec newCharges, int newNElectrons)
{
    basisHandler = newBasisHandler;
    integrator = new Integrator(basisHandler->getAngMomMax());
    nucleiPositions = newNucleiPositions;
    numberOfNuclei = nucleiPositions.n_rows;
    charges = newCharges;
    nElectrons = newNElectrons;
}



//--------------------------------------------------------------------------------------
// Returns one-electron integrals.
// 1st component = overlap
// 2nd component = kinetic energy + potential energy
rowvec2 System::getOneElectronIntegrals(int p, int q)
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
    int angMomTest = basisHandler->getAngMom(q);
    angMom = min(angMom, angMomTest);

    int nPrimitivesA = expA.n_elem;
    int nPrimitivesB = expB.n_elem;


    integrator->setPositionA(positionA);
    integrator->setPositionB(positionB);

    integrator->setMaxAngMom(angMom);

    int i = powersA(0);
    int j = powersB(0);
    int k = powersA(1);
    int l = powersB(1);
    int m = powersA(2);
    int n = powersB(2);

    double overlap = 0;
    double energy = 0;

    for (int v = 0; v < nPrimitivesA ; v++){
        for (int w = 0; w < nPrimitivesB; w++){

            integrator->setAlpha(expA(v));
            integrator->setBeta(expB(w));
            integrator->setE_AB();
            energy += integrator->kinetic(i, j, k, l, m, n);

            for (int x = 0; x < numberOfNuclei; x++){

                integrator->setPositionC(nucleiPositions.row(x));
                energy += -integrator->coulomb1(i,j,k,l,m,n)*charges(x);
            }

            energy *= coeffsA(v)*coeffsB(w);

            overlap += integrator->overlap(i,j,k,l,m,n)*coeffsA(v)*coeffsB(w);
        }
    }

    rowvec2 oneElectronIntegrals = {overlap, energy};
    return oneElectronIntegrals;
}



//------------------------------------------------------------------------------
// Returns two-electron integral
double System::getTwoElectronIntegral(int p, int r, int q, int s)
{
    rowvec expA = basisHandler->getExponents(p);
    rowvec expB = basisHandler->getExponents(q);
    rowvec expC = basisHandler->getExponents(r);
    rowvec expD = basisHandler->getExponents(s);

    rowvec coeffsA = basisHandler->getCoeffs(p);
    rowvec coeffsB = basisHandler->getCoeffs(q);
    rowvec coeffsC = basisHandler->getCoeffs(r);
    rowvec coeffsD = basisHandler->getCoeffs(s);

    irowvec powersA = basisHandler->getPowers(p);
    irowvec powersB = basisHandler->getPowers(q);
    irowvec powersC = basisHandler->getPowers(r);
    irowvec powersD = basisHandler->getPowers(s);

    rowvec3 positionA = basisHandler->getPosition(p);
    rowvec3 positionB = basisHandler->getPosition(q);
    rowvec3 positionC = basisHandler->getPosition(r);
    rowvec3 positionD = basisHandler->getPosition(s);

    irowvec indices = {p, r, q, s};
    int angMom = basisHandler->getAngMom(indices(0));
    int angMomTest;
    for (int i = 1; i < 4; i ++){
        angMomTest = basisHandler->getAngMom(indices(i));
        if (angMom < angMomTest){
            angMom = angMomTest;
        }
    }

    int nPrimitivesA = expA.n_elem;
    int nPrimitivesB = expB.n_elem;
    int nPrimitivesC = expC.n_elem;
    int nPrimitivesD = expD.n_elem;


    integrator->setPositionA(positionA);
    integrator->setPositionB(positionB);
    integrator->setPositionC(positionC);
    integrator->setPositionD(positionD);

    integrator->setMaxAngMom(angMom);

    int i1 = powersA(0);
    int j1 = powersB(0);
    int k1 = powersA(1);
    int l1 = powersB(1);
    int m1 = powersA(2);
    int n1 = powersB(2);
    int i2 = powersC(0);
    int j2 = powersD(0);
    int k2 = powersC(1);
    int l2 = powersD(1);
    int m2 = powersC(2);
    int n2 = powersD(2);

    double value = 0;

    for (int v = 0; v < nPrimitivesA; v++){
        for (int w = 0; w < nPrimitivesB; w++){
            for (int x = 0; x < nPrimitivesC; x++){
                for (int y = 0; y < nPrimitivesD; y++){

                    integrator->setAlpha(expA(v));
                    integrator->setBeta(expB(w));
                    integrator->setGamma(expC(x));
                    integrator->setDelta(expD(y));
                    integrator->setE_AB();
                    integrator->setE_CD();
                    value += integrator->coulomb2(i1,j1,k1,l1,m1,n1,i2,j2,k2,l2,m2,n2)
                            *coeffsA(v)*coeffsB(w)*coeffsC(x)*coeffsD(y);
                }
            }
        }
    }

    return value;
}


double System::getNucleiPotential()
{
    double value = 0;
    rowvec3 AB;

    for (int i = 0; i < numberOfNuclei; i++){
        for (int j = i+1; j < numberOfNuclei; j++){
            AB = nucleiPositions.row(i) - nucleiPositions.row(j);
            value += charges(i)*charges(j)/sqrt(dot(AB,AB));
        }
    }

    cout << value << endl;

    return value;
}


int System::getTotalNumOfBasisFunc()
{
    return basisHandler->getTotalNumOfBasisFunc();
}


int System::getNumOfElectrons()
{
    return nElectrons;
}

