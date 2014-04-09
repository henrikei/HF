#include "boysfunction.h"



BoysFunction::BoysFunction(int angMomMax)
{
    m_nMax = 4*angMomMax;
    m_F = zeros<vec>(4*angMomMax+1);
    m_Ftabulated.load("../../HartreeFock/inFiles/boys_tabulated.dat");
}


void BoysFunction::setx(double x)
{
    if (x <= 50){
        m_F(m_nMax) = tabulated(m_nMax, x);
    } else {
        m_F(m_nMax) = asymptotic(m_nMax, x);
    }

    double ex = exp(-x);

    for(int n = m_nMax; n > 0; n--){
        m_F(n-1) = (2*x*m_F(n) + ex)/(2*n - 1);
    }
}



double BoysFunction::returnValue(int n)
{
    return m_F(n);
}



double BoysFunction::tabulated(int n, double x)
{
    int nxVals = m_Ftabulated.n_rows;
    double dx = 50.0/(nxVals - 1); // x-spacing in tabulated values
    int xIndex = int ((x + 0.5*dx)/dx);
    double xt = xIndex*dx;     // tabulated x-value
    double Dx = x - xt;        // difference between actual and tabulated x-value

    double factorial = 1;
    double value = m_Ftabulated(xIndex, n);

    for(int k = 1; k < 7; k++){
        factorial *= k;
        value += m_Ftabulated(xIndex, n+k)*pow(-Dx,k)/factorial;
    }

    return value;
}



double BoysFunction::asymptotic(int n, double x)
{
    return factorial2(n)*sqrt(M_PI/(pow(x,2*n+1)))/(pow(2,n+1));
}

double BoysFunction::factorial2(int n)
{
    double value = 1;
    double i = 1;

    while(i < 2*n - 1){
        i += 2;
        value *= i;
    }

    return value;
}
