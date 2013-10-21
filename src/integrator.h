#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <armadillo>

using namespace arma;



class Integrator
{
public:
    Integrator();
    ~Integrator();
    void setAlpha(double newAlpha);
    void setBeta(double newBeta);
    void setPositionA(rowvec3 RA);
    void setPositionB(rowvec3 RB);
    void setMaxAngMom(int newAngMom);
    void setE();
    double overlap(int i, int j, int k, int l, int m, int n);
    double kinetic(int i, int j, int k, int l, int m, int n);

    cube E[3];

private:
    mat R;
    double alpha, beta;
    double angMom;
};

#endif // INTEGRATOR_H
