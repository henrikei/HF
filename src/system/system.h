#ifndef SYSTEM_H
#define SYSTEM_H

#include <armadillo>
#include "integrator/integrator.h"
#include "basisfunctions/basisfunctions.h"
#include "basisfunctions/contracted.h"
#include "basisfunctions/primitive.h"

using namespace std;
using namespace arma;



class System
{
public:
    System(BasisFunctions *basisFunctions, mat nucleiPositions, rowvec charges, int nElectrons);
    System(BasisFunctions *basisFunctions, mat nucleiPositions, rowvec charges, int nElectronsUp, int nElectronsDown);
    virtual ~System();

    rowvec2 getOneElectronIntegrals(int p, int q);
    double getTwoElectronIntegral(int p, int q, int r, int s);
    double getNucleiPotential();
    int getTotalNumOfBasisFunc();
    int getNumOfElectrons();
    int getNumOfElectronsUp();
    int getNumOfElectronsDown();
    mat getNucleiPositions();
    BasisFunctions* getBasisFunctions();

    void setNucleiPositions(mat nucleiPositions);

private:
    int m_nNuclei;
    mat m_nucleiPositions;
    rowvec m_charges;
    int m_nElectrons, m_nElectronsUp, m_nElectronsDown;
    BasisFunctions* m_basisFunctions;
    Integrator* m_integrator;
};

#endif // SYSTEM_H
