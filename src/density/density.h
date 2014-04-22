#ifndef DENSITY_H
#define DENSITY_H

#include <armadillo>
#include "basisfunctions/basisfunctions.h"

using namespace arma;

class Density
{
public:
    Density(BasisFunctions *basisFunctions, rowvec2 xlim, double dx,
            rowvec2 ylim, double dy, rowvec2 zlim, double dz);
    void setLimits(rowvec2 xlim, rowvec2 ylim, rowvec2 zlim);
    void setSpacing(double dx, double dy, double dz);
    virtual void printDensity(string filename)=0;

private:
    BasisFunctions *m_basisFunctions;
    rowvec2 m_xlim, m_ylim, m_zlim;
    double m_dx, m_dy, m_dz;
    cube m_density;

    void calcDensity(mat C, double factor);
};

#endif // DENSITY_H
