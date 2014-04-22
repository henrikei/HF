#ifndef CONTRACTED_H
#define CONTRACTED_H

#include <armadillo>
#include "primitive.h"

using namespace arma;
using namespace std;


class Contracted
{
public:
    Contracted(vector<Primitive*> primitives);
    int getNumOfPrimitives();
    Primitive *getPrimitive(int p);
    int getAngMom();
    double getValue(rowvec3 R);

private:
    vector<Primitive*> m_primitives;
};

#endif // CONTRACTED_H
