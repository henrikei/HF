#include "system.h"

System::System()
{

}

void System::setupBasis()
{
    double alpha = 20.2;
    double coeff = 1.0;
    char type = '1';

    PrimitiveBasis *primitive;
    primitive = new PrimitiveBasis(alpha, coeff, type);
    ContractedBasis *s1 = new ContractedBasis;
    s1->addPrimitive(primitive);
    basis.push_back(s1);
}
