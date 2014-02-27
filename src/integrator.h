#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <armadillo>
#include <vector>
#include "boysfunction.h"
#include "basisfunctions/primitive.h"

using namespace arma;
using namespace std;



class Integrator
{
public:
    Integrator(int angMomMax);
    ~Integrator();

    void setNucleusPosition(rowvec3 RC);
    void setPrimitiveA(Primitive* primitiveA);
    void setPrimitiveB(Primitive* primitiveB);
    void setPrimitiveC(Primitive* primitiveC);
    void setPrimitiveD(Primitive* primitiveD);
    void setE_AB(string intType);
    void setE_CD(string intType);
    double overlap();
    double kinetic();
    double coulomb1();
    double coulomb2();

private:
    Primitive *m_primitiveA, *m_primitiveB, *m_primitiveC, *m_primitiveD;
    cube E_AB[3];
    cube E_CD[3];
    vector<cube> R;
    BoysFunction *boys;
    rowvec3 nucleusPosition;
    double angMom;

    void setE(cube E[], Primitive *primitive1, Primitive *primitive2, string intType);
    void setR(double a, rowvec3 A, int tMax, int uMax, int vMax);
};

#endif // INTEGRATOR_H
