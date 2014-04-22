#include "contracted.h"

Contracted::Contracted(vector<Primitive *> primitives)
{
    m_primitives = primitives;
}

int Contracted::getNumOfPrimitives()
{
    return m_primitives.size();
}

Primitive *Contracted::getPrimitive(int p)
{
    return m_primitives.at(p);
}

int Contracted::getAngMom()
{
    return sum(m_primitives.at(0)->getPow());
}

double Contracted::getValue(rowvec3 R)
{
    rowvec3 center = m_primitives.at(0)->getPos();
    double xA = R(0) - center(0);
    double yA = R(1) - center(1);
    double zA = R(2) - center(2);
    double dist2 = xA*xA + yA*yA + zA*zA;

    int nPrim = m_primitives.size();
    double a, d;
    int i, j, k;
    double value = 0;

    for (int p = 0; p < nPrim; p++){
        a = m_primitives.at(p)->getExp();
        d = m_primitives.at(p)->getCoeff();
        i = m_primitives.at(p)->getPow()(0);
        j = m_primitives.at(p)->getPow()(1);
        k = m_primitives.at(p)->getPow()(2);

        value += d*pow(xA,i)*pow(yA,j)*pow(zA,k)*exp(-a*dist2);
    }
    return value;
}
