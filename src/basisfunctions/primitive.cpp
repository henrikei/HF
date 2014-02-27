#include "primitive.h"

Primitive::Primitive(double exp, double coeff, irowvec pow, rowvec3 position)
{
    m_exp = exp;
    m_coeff = coeff;
    m_pow = pow;
    m_position = position;
}

double Primitive::getExp()
{
    return m_exp;
}

double Primitive::getCoeff()
{
    return m_coeff;
}

irowvec3 Primitive::getPow()
{
    return m_pow;
}

rowvec3 Primitive::getPos()
{
    return m_position;
}
