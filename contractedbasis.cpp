#include "contractedbasis.h"

ContractedBasis::ContractedBasis()
{
}

void ContractedBasis::addPrimitive(PrimitiveBasis *primitive)
{
    m_primitives.push_back(primitive);
}

vector<PrimitiveBasis*> ContractedBasis::getPrimitives()
{
    return m_primitives;
}
