#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#include <armadillo>

using namespace arma;


class BoysFunction
{
public:
    BoysFunction(int angMomMax);
    void setx(double x);
    double returnValue(int n);
private:
    void setFs(double x);
    double tabulated(int n, double x);
    double asymptotic(int n, double x);
    double factorial2(int n);
    void readFromFile();

    mat m_Ftabulated;
    vec m_F;
    int m_nMax;
    int m_nx;
    int m_nN;
};

#endif // BOYSFUNCTION_H
