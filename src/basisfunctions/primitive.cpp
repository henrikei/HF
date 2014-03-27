#include "primitive.h"

Primitive::Primitive(double exp, double coeff, irowvec pow, int posNum)
{
    m_exp = exp;
    m_coeff = coeff;
    m_pow = pow;
    m_posNum = posNum;
    m_nucleiPositions = NULL;
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
    return m_nucleiPositions->row(m_posNum);
}

void Primitive::setPosPointer(mat *nucleiPositions)
{
    m_nucleiPositions = nucleiPositions;
}
