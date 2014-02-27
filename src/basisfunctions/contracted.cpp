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
