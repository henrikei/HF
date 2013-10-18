#include "primitivebasis.h"

PrimitiveBasis::PrimitiveBasis(double alpha, double d, char type)
{
    m_alpha = alpha;
    m_coeff = d;
    m_type = type;
}

double PrimitiveBasis::getAlpha()
{
    return m_alpha;
}

double PrimitiveBasis::getCoeff()
{
    return m_coeff;
}
