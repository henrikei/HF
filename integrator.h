#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <armadillo>

using namespace arma;



class Integrator
{
public:
    Integrator();
    ~Integrator();
    void setE();
    mat R;
    double a, b;
    double l;

    cube E[3];
};

#endif // INTEGRATOR_H
