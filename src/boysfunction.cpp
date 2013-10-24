#include "boysfunction.h"



BoysFunction::BoysFunction(int nMaxNew, double x)
{
    nMax = nMaxNew;
    F = zeros<vec>(nMax+1);
    readFromFile();
    setFs(x);
}

void BoysFunction::resetx(double x)
{
    setFs(x);
}



double BoysFunction::returnValue(int n)
{
    return F(n);
}



void BoysFunction::setFs(double x)
{
    if (x <= 30){
        F(nMax) = tabulated(nMax, x);
    } else {
        F(nMax) = asymptotic(nMax, x);
    }

    double ex = exp(-x);

    for(int n = nMax; n > 0; n--){
        F(n-1) = (2*x*F(n) + ex)/(2*n - 1);
    }
}



double BoysFunction::tabulated(int n, double x)
{
    double dx = 30.0/(nx - 1); // x-spacing in tabulated values
    int xIndex = int ((x + 0.5*dx)/dx);
    double xt = xIndex*dx;     // tabulated x-value
    double Dx = x - xt;        // difference between actual and tabulated x-value

    double value = 0;
    double factorial = 1;

    for(int k = 0; k < 6; k++){
        if(k != 0){
            factorial *= k;
        }
        value += Ftabulated(xIndex, n+k)*pow(-Dx,k)/factorial;
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



void BoysFunction::readFromFile()
{
    Ftabulated.load("boys_tabulated.dat");
    nx = Ftabulated.n_rows;
    nN = Ftabulated.n_cols;
}
