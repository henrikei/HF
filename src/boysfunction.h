#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#include <armadillo>

using namespace arma;


class BoysFunction
{
public:
    BoysFunction(int nMaxNew, double x);
    void resetx(double x);
    double returnValue(int n);
private:
    void setFs(double x);
    double tabulated(int n, double x);
    double asymptotic(int n, double x);
    double factorial2(int n);
    void readFromFile();

    mat Ftabulated;
    vec F;
    int nMax;
    int nx;
    int nN;
};

#endif // BOYSFUNCTION_H
