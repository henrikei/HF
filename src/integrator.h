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
    virtual ~Integrator();

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

    // The following are public only for unit testing
    double overlap(int i, int j, int k, int l, int m, int n);
    double kinetic(int i, int j, int k, int l, int m, int n);
    double coulomb1(int i, int j, int k, int l, int m, int n);

private:
    Primitive *m_primitiveA, *m_primitiveB, *m_primitiveC, *m_primitiveD;
    cube m_E_AB[3];
    cube m_E_CD[3];
    vector<cube> m_R;
    BoysFunction *m_boys;
    rowvec3 m_nucleusPosition;
    double m_angMom;

    void setE(cube E[], Primitive *primitive1, Primitive *primitive2, string intType);
    void setR(double a, rowvec3 A, int tMax, int uMax, int vMax);
};

#endif // INTEGRATOR_H
