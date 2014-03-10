#ifndef TWODIMTEST_H
#define TWODIMTEST_H

#include <armadillo>
#include "func.h"

using namespace std;
using namespace arma;


class TwoDimTest : public Func
{
public:
    TwoDimTest();
    virtual double getValue(rowvec x);
};

#endif // TWODIMTEST_H
