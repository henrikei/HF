#include "system.h"



System::System(BasisFunctions *basisFunctions, mat nucleiPositions, rowvec charges, int nElectrons)
{
    m_basisFunctions = basisFunctions;
    m_integrator = new Integrator(basisFunctions->getAngMomMax());
    m_nucleiPositions = nucleiPositions;
    m_charges = charges;
    m_nElectrons = nElectrons;
    m_nElectronsUp = nElectrons/2;
    m_nElectronsDown = nElectrons/2;
    // This constructor places +1 electron spin-up if number of electrons is odd
    if (m_nElectrons % 2 == 1){
        m_nElectronsUp += 1;
    }

    m_basisFunctions->setPosPointer(&m_nucleiPositions);
    m_nNuclei = m_nucleiPositions.n_rows;
}

System::System(BasisFunctions *basisFunctions, mat nucleiPositions, rowvec charges, int nElectronsUp, int nElectronsDown):
    System(basisFunctions, nucleiPositions, charges, nElectronsUp + nElectronsDown)
{
    m_nElectronsUp = nElectronsUp;
    m_nElectronsDown = nElectronsDown;
}

System::~System()
{
    delete m_integrator;
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

    Primitive* primitiveA;
    Primitive* primitiveB;
    for (int v = 0; v < nPrimitivesA ; v++){
        primitiveA = contractedA->getPrimitive(v);
        m_integrator->setPrimitiveA(primitiveA);
        for (int w = 0; w < nPrimitivesB; w++){
            primitiveB = contractedB->getPrimitive(w);
            m_integrator->setPrimitiveB(primitiveB);

            m_integrator->setE_AB("oneParticle");

            energy += m_integrator->kinetic()*primitiveA->getCoeff()*primitiveB->getCoeff();

            for (int x = 0; x < m_nNuclei; x++){

                m_integrator->setNucleusPosition(m_nucleiPositions.row(x));
                energy += -m_integrator->coulomb1()*m_charges(x)*primitiveA->getCoeff()*primitiveB->getCoeff();
            }

            overlap += m_integrator->overlap()*primitiveA->getCoeff()*primitiveB->getCoeff();
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

    Primitive* primitiveA;
    Primitive* primitiveB;
    Primitive* primitiveC;
    Primitive* primitiveD;

    for (int v = 0; v < nPrimitivesA; v++){

        primitiveA = contractedA->getPrimitive(v);
        m_integrator->setPrimitiveA(primitiveA);

        for (int w = 0; w < nPrimitivesB; w++){

            primitiveB = contractedB->getPrimitive(w);
            m_integrator->setPrimitiveB(primitiveB);
            m_integrator->setE_AB("twoParticle");

            for (int x = 0; x < nPrimitivesC; x++){

                primitiveC = contractedC->getPrimitive(x);
                m_integrator->setPrimitiveC(primitiveC);

                for (int y = 0; y < nPrimitivesD; y++){

                    primitiveD = contractedD->getPrimitive(y);
                    m_integrator->setPrimitiveD(primitiveD);
                    m_integrator->setE_CD("twoParticle");

                    value += m_integrator->coulomb2()
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

    for (int i = 0; i < m_nNuclei; i++){
        for (int j = i+1; j < m_nNuclei; j++){
            AB = m_nucleiPositions.row(i) - m_nucleiPositions.row(j);
            value += m_charges(i)*m_charges(j)/sqrt(dot(AB,AB));
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
    return m_nElectrons;
}

int System::getNumOfElectronsUp()
{
    return m_nElectronsUp;
}

int System::getNumOfElectronsDown()
{
    return m_nElectronsDown;
}

mat System::getNucleiPositions()
{
    return m_nucleiPositions;
}

BasisFunctions *System::getBasisFunctions()
{
    return m_basisFunctions;
}

void System::setNucleiPositions(mat nucleiPositions)
{
    m_nucleiPositions = nucleiPositions;
}

