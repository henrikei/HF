#include "twodimtest.h"

TwoDimTest::TwoDimTest(rowvec xInitial)
{
    m_xInitial = xInitial;
}

rowvec TwoDimTest::getxInitial()
{
    return m_xInitial;
}

double TwoDimTest::getValue(rowvec x)
{
    return 100*(x(1) - x(0)*x(0))*(x(1) - x(0)*x(0)) + (2.4 - x(0))*(2.4 - x(0));
}
