#include "integrator.h"

Integrator::Integrator()
{
    R = zeros(2,3);

    alpha = 0;
    beta = 0;
    angMom = 0;
}

Integrator::~Integrator()
{
}

void Integrator::setAlpha(double newAlpha)
{
    alpha = newAlpha;
}

void Integrator::setBeta(double newBeta)
{
    beta = newBeta;
}

void Integrator::setPositionA(rowvec3 RA)
{
    R.row(0) = RA;
}

void Integrator::setPositionB(rowvec3 RB)
{
    R.row(1) = RB;
}

void Integrator::setMaxAngMom(int newAngMom)
{
    angMom = newAngMom;
}

void Integrator::setE()
{
    double p = alpha + beta;
    rowvec3 P = (alpha*R.row(0) + beta*R.row(1))/p;
    rowvec3 AB = R.row(0) - R.row(1);
    rowvec3 PA = P - R.row(0);
    rowvec3 PB = P - R.row(1);

    E[0] = zeros<cube>(angMom+3,angMom+3,2*(angMom+3)+1);
    E[1] = zeros<cube>(angMom+3,angMom+3,2*(angMom+3)+1);
    E[2] = zeros<cube>(angMom+3,angMom+3,2*(angMom+3)+1);


    // Loop over x-, y- and z- coordinates

    for(int dir = 0; dir < 3; dir++){

        // First loop over t and i with j = 0

        E[dir](0,0,0) = exp(-alpha*beta*AB(dir)*AB(dir)/p);
        for (int i = 0; i < angMom+2; i++){
            // Treat the case t=0 separately due to (t-1) term
            E[dir](i+1,0,0) = PA(dir)*E[dir](i,0,0) + E[dir](i,0,1);
            for (int t = 1; t <= i + 0 + 1; t++){
                E[dir](i+1,0,t) = E[dir](i,0,t-1)/(2*p) + PA(dir)*E[dir](i,0,t) + (t+1)*E[dir](i,0,t+1);
            }
        }

        // Second loop over t and j and i

        // Must here let i <= l because the forward loop is on index j
        for (int i = 0; i <= angMom+2; i++){
            for (int j = 0; j < angMom+2; j++){
                E[dir](i,j+1,0) = PB(dir)*E[dir](i,j,0) + E[dir](i,j,1);
                for (int t = 1; t <= i + j + 1; t++){
                    E[dir](i,j+1,t) = E[dir](i,j,t-1)/(2*p) + PB(dir)*E[dir](i,j,t) + (t+1)*E[dir](i,j,t+1);
                }
            }
        }
    }
}

double Integrator::overlap(int i, int j, int k, int l, int m, int n)
{
    double p = alpha + beta;
    return E[0](i,j,0)*E[1](k,l,0)*E[2](m,n,0)*pow(M_PI/p, 1.5);
}

double Integrator::kinetic(int i, int j, int k, int l, int m, int n)
{
    double Tij, Tkl, Tmn, kinetic;
    double p = alpha + beta;

    if (j < 2){
        Tij = 4*beta*beta*E[0](i,j+2,0) - 2*beta*(2*j + 1)*E[0](i,j,0);
    } else {
        Tij = 4*beta*beta*E[0](i,j+2,0) - 2*beta*(2*j + 1)*E[0](i,j,0) + j*(j-1)*E[0](i,j-2,0);
    }

    if (l < 2){
        Tkl = 4*beta*beta*E[1](k,l+2,0) - 2*beta*(2*l + 1)*E[1](k,l,0);
    } else {
        Tkl = 4*beta*beta*E[1](k,l+2,0) - 2*beta*(2*l + 1)*E[1](k,l,0) + l*(l-1)*E[1](k,l-2,0);
    }

    if (n < 2){
        Tmn = 4*beta*beta*E[2](m,n+2,0) - 2*beta*(2*n + 1)*E[2](m,n,0);
    } else {
        Tmn = 4*beta*beta*E[2](m,n+2,0) - 2*beta*(2*n + 1)*E[2](m,n,0) + n*(n-1)*E[2](m,n-2,0);
    }

    kinetic = -0.5*(Tij*E[1](k,l,0)*E[2](m,n,0) + E[0](i,j,0)*Tkl*E[2](m,n,0) + E[0](i,j,0)*E[1](k,l,0)*Tmn)*pow(M_PI/p, 1.5);

    return kinetic;
}
