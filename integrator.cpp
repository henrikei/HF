#include "integrator.h"

Integrator::Integrator()
{
    R = zeros(2,3);
    rowvec RA = {1.2, 2.3, 3.4};
    R.row(0) = RA;
    rowvec RB = {-1.3, 1.4, -2.4};
    R.row(1) = RB;

    a = 0.2;
    b = 0.3;
    l = 2;

    E[0] = zeros<cube>(l+1,l+1,2*(l+1)+1);
    E[1] = zeros<cube>(l+1,l+1,2*(l+1)+1);
    E[2] = zeros<cube>(l+1,l+1,2*(l+1)+1);
}

Integrator::~Integrator()
{
}

void Integrator::setE()
{
    double p = a + b;
    rowvec3 P = (a*R.row(0) + b*R.row(1))/p;
    rowvec3 AB = R.row(0) - R.row(1);
    rowvec3 PA = P - R.row(0);
    rowvec3 PB = P - R.row(1);


    // Loop over x-, y- and z- coordinates

    for(int dir = 0; dir < 3; dir++){

        // First loop over t and i with j = 0

        E[dir](0,0,0) = exp(-a*b*AB(dir)*AB(dir)/p);
        for (int i = 0; i < l; i++){
            // Treat the case t=0 separately due to (t-1) term
            E[dir](i+1,0,0) = PA(dir)*E[dir](i,0,0) + E[dir](i,0,1);
            for (int t = 1; t <= i + 0 + 1; t++){
                E[dir](i+1,0,t) = E[dir](i,0,t-1)/(2*p) + PA(dir)*E[dir](i,0,t) + (t+1)*E[dir](i,0,t+1);
            }
        }

        // Second loop over t and j and i

        // Must here let i <= l because the forward loop is on index j
        for (int i = 0; i <= l; i++){
            for (int j = 0; j < l; j++){
                E[dir](i,j+1,0) = PB(dir)*E[dir](i,j,0) + E[dir](i,j,1);
                for (int t = 1; t <= i + j + 1; t++){
                    E[dir](i,j+1,t) = E[dir](i,j,t-1)/(2*p) + PB(dir)*E[dir](i,j,t) + (t+1)*E[dir](i,j,t+1);
                }
            }
        }
    }
}
