#include "density.h"

Density::Density(BasisFunctions *basisFunctions, rowvec2 xlim, double dx, rowvec2 ylim, double dy, rowvec2 zlim, double dz)
{
    m_basisFunctions = basisFunctions;
    m_xlim = xlim; m_ylim = ylim; m_zlim = zlim;
    m_dx = dx; m_dy = dy; m_dz = dz;

    int nx = (m_xlim(1) - m_xlim(0))/m_dx + 1;
    int ny = (m_ylim(1) - m_ylim(0))/m_dy + 1;
    int nz = (m_zlim(1) - m_zlim(0))/m_dz + 1;
    m_density = zeros(nx, ny, nz);
}

void Density::setLimits(rowvec2 xlim, rowvec2 ylim, rowvec2 zlim)
{
    m_xlim = xlim; m_ylim = ylim; m_zlim = zlim;
}

void Density::setSpacing(double dx, double dy, double dz)
{
    m_dx = dx; m_dy = dy; m_dz = dz;
}

void Density::calcDensity(mat C, double factor)
{
}
