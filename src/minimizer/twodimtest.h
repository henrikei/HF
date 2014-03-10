#ifndef TWODIMTEST_H
#define TWODIMTEST_H

#include <armadillo>
#include "func.h"

using namespace std;
using namespace arma;


class TwoDimTest : public Func
{
public:
    TwoDimTest(rowvec xInitial);
    virtual rowvec getxInitial();
    virtual double getValue(rowvec x);
private:
    rowvec m_xInitial;
};

#endif // TWODIMTEST_H
