#include "integrator.h"

Integrator::Integrator()
{
    Rnuclei = zeros(2,3);

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
    Rnuclei.row(0) = RA;
}

void Integrator::setPositionB(rowvec3 RB)
{
    Rnuclei.row(1) = RB;
}

void Integrator::setMaxAngMom(int newAngMom)
{
    angMom = newAngMom;
}

void Integrator::setE()
{
    double p = alpha + beta;
    rowvec3 P = (alpha*Rnuclei.row(0) + beta*Rnuclei.row(1))/p;
    rowvec3 AB = Rnuclei.row(0) - Rnuclei.row(1);
    rowvec3 PA = P - Rnuclei.row(0);
    rowvec3 PB = P - Rnuclei.row(1);

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

void Integrator::setR(rowvec3 C)
{
    double p = alpha + beta;
    rowvec3 P = (alpha*Rnuclei.row(0) + beta*Rnuclei.row(1))/p;
    rowvec3 PC = P - C;

    int tMax = 2*angMom;
    int nMax = 3*tMax;
    BoysFunction Boys(nMax, p*dot(PC,PC));

    cube Rinit = zeros<cube>(tMax+1, tMax+1, tMax+1);

    // Generate initial Rs
    for(int n = 0; n < nMax+1; n++){
        Rinit(0,0,0) = pow(-2*p, n)*Boys.returnValue(n);
        R.push_back(Rinit);
    }


    for(int tuv = 1; tuv < nMax+1; tuv++){
        for(int n = 0; n < tuv+1; n++){
            for(int t = 0; t < tMax; t++){
                for(int u = 0; u < tMax; u++){
                    for(int v = 0; v < tMax; v++){
                        if(t + u + v == tuv){

                            if(t == 0){
                                R.at(n)(t+1,u,v) = PC(0)*R.at(n+1)(t,u,v);
                            } else {
                                R.at(n)(t+1,u,v) = t*R.at(n+1)(t-1,u,v) + PC(0)*R.at(n+1)(t,u,v);
                            }

                            if(u == 0){
                                R.at(n)(t,u+1,v) = PC(1)*R.at(n+1)(t,u,v);
                            } else {
                                R.at(n)(t,u+1,v) = u*R.at(n+1)(t,u-1,v) + PC(1)*R.at(n+1)(t,u,v);
                            }

                            if(v == 0){
                                R.at(n)(t,u,v+1) = PC(2)*R.at(n+1)(t,u,v);
                            } else {
                                R.at(n)(t,u,v+1) = v*R.at(n+1)(t,u,v-1) + PC(2)*R.at(n+1)(t,u,v);
                            }

                        }
                    }
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

double Integrator::coulomb1(int i, int j, int k, int l, int m, int n)
{
    int tMax = i + j;
    int uMax = k + l;
    int vMax = m + n;
    double p = alpha + beta;

    double value = 0;

    for (int t = 0; t < tMax + 1; t++){
        for (int u = 0; u < uMax + 1; u++){
            for (int v = 0; v < vMax + 1; v++){
                value += E[0](i,j,t)*E[1](k,l,u)*E[2](m,n,v)*R.at(0)(t,u,v);
            }
        }
    }

    value *= 2*M_PI/p;

    return value;
}
