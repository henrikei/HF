#ifndef FUNC_H
#define FUNC_H

#include <armadillo>

using namespace std;
using namespace arma;


class Func
{
public:
    Func();
    virtual double getValue(rowvec x)=0;
};

#endif // FUNC_H
