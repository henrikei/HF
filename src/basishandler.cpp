#include "basishandler.h"

BasisHandler::BasisHandler()
{
}


void BasisHandler::addBasisFunctions(BasisFunctions *basis)
{
    int nucleiNumber = allBasisFunctions.size();
    allBasisFunctions.push_back(basis);
    int n = basis->getNumOfBasisFunc();

    for (int i = 0; i < n; i++){
        map.push_back(nucleiNumber);
    }
}


int BasisHandler::getTotalNumOfBasisFunc()
{
    return map.size();
}


rowvec3 BasisHandler::getPosition(int p)
{
    int nucleiNumber = map.at(p);
    return allBasisFunctions.at(nucleiNumber)->getPosition();
}


rowvec BasisHandler::getExponents(int p)
{
    int n = 0;
    int nucleiNumber = map.at(p);

    for (int i = 0; i < nucleiNumber; i++){
        n += allBasisFunctions.at(i)->getNumOfBasisFunc();
    }

    return allBasisFunctions.at(nucleiNumber)->getExponents(p-n);
}


rowvec BasisHandler::getCoeffs(int p)
{
    int n = 0;
    int nucleiNumber = map.at(p);

    for (int i = 0; i < nucleiNumber; i++){
        n += allBasisFunctions.at(i)->getNumOfBasisFunc();
    }

    return allBasisFunctions.at(nucleiNumber)->getCoeffs(p-n);
}


irowvec BasisHandler::getPowers(int p)
{
    int n = 0;
    int nucleiNumber = map.at(p);

    for (int i = 0; i < nucleiNumber; i++){
        n += allBasisFunctions.at(i)->getNumOfBasisFunc();
    }

    return allBasisFunctions.at(nucleiNumber)->getPowers(p-n);
}

int BasisHandler::getAngMom(int p)
{
    int nucleiNumber = map.at(p);

    return allBasisFunctions.at(nucleiNumber)->getAngMom();
}


int BasisHandler::getAngMomMax()
{
    int angMomMax = 0;
    int angMomTest = 0;
    int nNuclei = allBasisFunctions.size();

    for (int i = 0; i < nNuclei; i++){
        angMomTest = allBasisFunctions.at(i)->getAngMom();
        if (angMomTest > angMomMax){
            angMomMax = angMomTest;
        }
    }

    return angMomMax;
}
