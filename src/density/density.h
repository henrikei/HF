#ifndef DENSITY_H
#define DENSITY_H

#include <armadillo>
#include <fstream>
#include "basisfunctions/basisfunctions.h"

using namespace arma;

class Density
{
public:
    Density(BasisFunctions *basisFunctions, field<mat> P, rowvec3 R1, rowvec3 R2,
            double dx, double dy, double dz);
    void printDensity(string filename);

private:
    BasisFunctions *m_basisFunctions;
    field<mat> m_P;
    rowvec3 m_R1, m_R2;
    double m_dx, m_dy, m_dz;

    void printDensity(mat P, string filename);
};

#endif // DENSITY_H
