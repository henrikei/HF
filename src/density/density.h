#ifndef DENSITY_H
#define DENSITY_H

#include <armadillo>
#include <fstream>
#include "basisfunctions/basisfunctions.h"

using namespace arma;

class Density
{
public:
    Density(BasisFunctions *basisFunctions, rowvec3 R1, rowvec3 R2,
            double dx, double dy, double dz);
    void printDensity(field<mat> P, string filename);
    void printSpinDensity(field<mat> P, string filename);
    void printMolecularOrbital(field<mat> C, int orbital, string filename, int spin=0);

private:
    BasisFunctions *m_basisFunctions;
    rowvec3 m_R1, m_R2;
    double m_dx, m_dy, m_dz;

    void printDensity(mat P, string filename);
};

#endif // DENSITY_H
