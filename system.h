#ifndef SYSTEM_H
#define SYSTEM_H

#include "contractedbasis.h"

class System
{
public:
    System();
    void setupBasis();
private:
    vector<ContractedBasis*> basis;
};

#endif // SYSTEM_H
