#include "twodimtest.h"

TwoDimTest::TwoDimTest()
{
}

double TwoDimTest::getValue(rowvec x)
{
    return x(0)*x(0) + x(1)*x(1) + 2*x(0) - 4*x(1) + 10;
}
