#include "system.h"



System::System(BasisFunctions2* basisFunctions, mat newNucleiPositions, rowvec newCharges, int newNElectrons)
{
    m_basisFunctions = basisFunctions;
    integrator = new Integrator(basisFunctions->getAngMomMax());
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
    Contracted* contractedA = m_basisFunctions->getContracted(p);
    Contracted* contractedB = m_basisFunctions->getContracted(q);

    int nPrimitivesA = contractedA->getNumOfPrimitives();
    int nPrimitivesB = contractedB->getNumOfPrimitives();

    double overlap = 0;
    double energy = 0;

    for (int v = 0; v < nPrimitivesA ; v++){
        Primitive* primitiveA = contractedA->getPrimitive(v);
        for (int w = 0; w < nPrimitivesB; w++){
            Primitive* primitiveB = contractedB->getPrimitive(w);

            integrator->setPrimitiveA(primitiveA);
            integrator->setPrimitiveB(primitiveB);

            integrator->setE_AB("oneParticle");

            energy += integrator->kinetic()*primitiveA->getCoeff()*primitiveB->getCoeff();

            for (int x = 0; x < numberOfNuclei; x++){

                integrator->setNucleusPosition(nucleiPositions.row(x));
                energy += -integrator->coulomb1()*charges(x)*primitiveA->getCoeff()*primitiveB->getCoeff();
            }

            overlap += integrator->overlap()*primitiveA->getCoeff()*primitiveB->getCoeff();
        }
    }

    rowvec2 oneElectronIntegrals = {overlap, energy};
    return oneElectronIntegrals;
}



//------------------------------------------------------------------------------
// Returns two-electron integral
double System::getTwoElectronIntegral(int p, int q, int r, int s)
{
    Contracted* contractedA = m_basisFunctions->getContracted(p);
    Contracted* contractedB = m_basisFunctions->getContracted(r);
    Contracted* contractedC = m_basisFunctions->getContracted(q);
    Contracted* contractedD = m_basisFunctions->getContracted(s);

    int nPrimitivesA = contractedA->getNumOfPrimitives();
    int nPrimitivesB = contractedB->getNumOfPrimitives();
    int nPrimitivesC = contractedC->getNumOfPrimitives();
    int nPrimitivesD = contractedD->getNumOfPrimitives();

    double value = 0;

    for (int v = 0; v < nPrimitivesA; v++){
        Primitive* primitiveA = contractedA->getPrimitive(v);
        for (int w = 0; w < nPrimitivesB; w++){
            Primitive* primitiveB = contractedB->getPrimitive(w);
            for (int x = 0; x < nPrimitivesC; x++){
                Primitive* primitiveC = contractedC->getPrimitive(x);
                for (int y = 0; y < nPrimitivesD; y++){
                    Primitive* primitiveD = contractedD->getPrimitive(y);

                    integrator->setPrimitiveA(primitiveA);
                    integrator->setPrimitiveB(primitiveB);
                    integrator->setPrimitiveC(primitiveC);
                    integrator->setPrimitiveD(primitiveD);

                    integrator->setE_AB("twoParticle");
                    integrator->setE_CD("twoParticle");
                    value += integrator->coulomb2()
                            *primitiveA->getCoeff()*primitiveB->getCoeff()*primitiveC->getCoeff()*primitiveD->getCoeff();
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

    return value;
}


int System::getTotalNumOfBasisFunc()
{
    return m_basisFunctions->getNumOfContracteds();
}


int System::getNumOfElectrons()
{
    return nElectrons;
}

