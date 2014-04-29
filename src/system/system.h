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
    System(BasisFunctions* m_basisFunctions, mat nucleiPositions, rowvec charges, int nElectrons);
    virtual ~System();

    rowvec2 getOneElectronIntegrals(int p, int q);
    double getTwoElectronIntegral(int p, int q, int r, int s);
    double getNucleiPotential();
    int getTotalNumOfBasisFunc();
    int getNumOfElectrons();
    mat getNucleiPositions();
    BasisFunctions* getBasisFunctions();

    void setNucleiPositions(mat nucleiPositions);

private:
    int m_nNuclei;
    mat m_nucleiPositions;
    rowvec m_charges;
    int m_nElectrons;
    BasisFunctions* m_basisFunctions;
    Integrator* m_integrator;
};

#endif // SYSTEM_H
