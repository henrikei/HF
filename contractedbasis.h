#ifndef CONTRACTEDBASIS_H
#define CONTRACTEDBASIS_H

#include "primitivebasis.h"
#include <vector>

using namespace std;

class ContractedBasis
{
public:
    ContractedBasis();
    void addPrimitive(PrimitiveBasis *primitive);
    vector <PrimitiveBasis*> getPrimitives();
private:
    vector<PrimitiveBasis*> m_primitives;
};

#endif // CONTRACTEDBASIS_H
