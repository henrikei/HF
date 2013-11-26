#include "integrator.h"

Integrator::Integrator(int angMomMax)
{
    boys = new BoysFunction(angMomMax);

    Rnuclei = zeros(4,3);

    alpha = 0;
    beta = 0;
    gamma = 0;
    delta = 0;
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

void Integrator::setGamma(double newGamma)
{
    gamma = newGamma;
}

void Integrator::setDelta(double newDelta)
{
    delta = newDelta;
}

void Integrator::setPositionA(rowvec3 RA)
{
    Rnuclei.row(0) = RA;
}

void Integrator::setPositionB(rowvec3 RB)
{
    Rnuclei.row(1) = RB;
}

void Integrator::setPositionC(rowvec3 RC)
{
    Rnuclei.row(2) = RC;
}

void Integrator::setPositionD(rowvec3 RD)
{
    Rnuclei.row(3) = RD;
}

void Integrator::setMaxAngMom(int newAngMom)
{
    angMom = newAngMom;
}

void Integrator::setE_AB()
{
    double p = alpha + beta;
    rowvec3 P = (alpha*Rnuclei.row(0) + beta*Rnuclei.row(1))/p;
    rowvec3 AB = Rnuclei.row(0) - Rnuclei.row(1);
    rowvec3 PA = P - Rnuclei.row(0);
    rowvec3 PB = P - Rnuclei.row(1);

    E_AB[0] = zeros<cube>(angMom+3,angMom+3,2*(angMom+2)+2);
    E_AB[1] = zeros<cube>(angMom+3,angMom+3,2*(angMom+2)+2);
    E_AB[2] = zeros<cube>(angMom+3,angMom+3,2*(angMom+2)+2);


    // Loop over x-, y- and z- coordinates

    for(int dir = 0; dir < 3; dir++){

        // First loop over t and i with j = 0

        E_AB[dir](0,0,0) = exp(-alpha*beta*AB(dir)*AB(dir)/p);
        for (int i = 0; i < angMom+2; i++){
            // Treat the case t=0 separately due to (t-1) term
            E_AB[dir](i+1,0,0) = PA(dir)*E_AB[dir](i,0,0) + E_AB[dir](i,0,1);
            for (int t = 1; t <= i + 0 + 1; t++){
                E_AB[dir](i+1,0,t) = E_AB[dir](i,0,t-1)/(2*p) + PA(dir)*E_AB[dir](i,0,t) + (t+1)*E_AB[dir](i,0,t+1);
            }
        }

        // Second loop over t and j and i

        // Must here let i <= l because the forward loop is on index j
        for (int i = 0; i <= angMom+2; i++){
            for (int j = 0; j < angMom+2; j++){
                E_AB[dir](i,j+1,0) = PB(dir)*E_AB[dir](i,j,0) + E_AB[dir](i,j,1);
                for (int t = 1; t <= i + j + 1; t++){
                    E_AB[dir](i,j+1,t) = E_AB[dir](i,j,t-1)/(2*p) + PB(dir)*E_AB[dir](i,j,t) + (t+1)*E_AB[dir](i,j,t+1);
                }
            }
        }
    }
}



void Integrator::setE_CD()
{
    double q = gamma + delta;
    rowvec3 Q = (gamma*Rnuclei.row(2) + delta*Rnuclei.row(3))/q;
    rowvec3 CD = Rnuclei.row(2) - Rnuclei.row(3);
    rowvec3 QC = Q - Rnuclei.row(2);
    rowvec3 QD = Q - Rnuclei.row(3);

    E_CD[0] = zeros<cube>(angMom+3,angMom+3,2*(angMom+2)+2);
    E_CD[1] = zeros<cube>(angMom+3,angMom+3,2*(angMom+2)+2);
    E_CD[2] = zeros<cube>(angMom+3,angMom+3,2*(angMom+2)+2);


    // Loop over x-, y- and z- coordinates

    for(int dir = 0; dir < 3; dir++){

        // First loop over t and i with j = 0

        E_CD[dir](0,0,0) = exp(-gamma*delta*CD(dir)*CD(dir)/q);
        for (int i = 0; i < angMom+2; i++){
            // Treat the case t=0 separately due to (t-1) term
            E_CD[dir](i+1,0,0) = QC(dir)*E_CD[dir](i,0,0) + E_CD[dir](i,0,1);
            for (int t = 1; t <= i + 0 + 1; t++){
                E_CD[dir](i+1,0,t) = E_CD[dir](i,0,t-1)/(2*q) + QC(dir)*E_CD[dir](i,0,t) + (t+1)*E_CD[dir](i,0,t+1);
            }
        }

        // Second loop over t and j and i

        // Must here let i <= l because the forward loop is on index j
        for (int i = 0; i <= angMom+2; i++){
            for (int j = 0; j < angMom+2; j++){
                E_CD[dir](i,j+1,0) = QD(dir)*E_CD[dir](i,j,0) + E_CD[dir](i,j,1);
                for (int t = 1; t <= i + j + 1; t++){
                    E_CD[dir](i,j+1,t) = E_CD[dir](i,j,t-1)/(2*q) + QD(dir)*E_CD[dir](i,j,t) + (t+1)*E_CD[dir](i,j,t+1);
                }
            }
        }
    }
}



void Integrator::setR(double a, rowvec3 A, int tMax, int uMax, int vMax)
{
//    int tMax, nMax;

//    if (coulombType == 0){
//        tMax = 2*angMom;
//        nMax = 2*angMom;
//    } else {
//        tMax = 4*angMom;
//        nMax = 4*angMom;
//    }

    int nMax = tMax + uMax + vMax;

    boys->setx(a*dot(A,A));

    cube Rinit = zeros<cube>(tMax+1, uMax+1, vMax+1);

    // Generate initial Rs
    R.clear();
    for(int n = 0; n < nMax+1; n++){
        Rinit(0,0,0) = pow(-2*a, n)*boys->returnValue(n);
        R.push_back(Rinit);
    }


    for(int t = 0; t < tMax; t++){
        for(int n = 0; n < nMax - t; n++){
            if (t == 0){
                R.at(n)(t+1,0,0) = A(0)*R.at(n+1)(t,0,0);
            } else {
                R.at(n)(t+1,0,0) = t*R.at(n+1)(t-1,0,0) + A(0)*R.at(n+1)(t,0,0);
            }
        }
    }

    for(int t = 0; t < tMax+1; t++){
        for(int u = 0; u < uMax; u++){
            for(int n = 0; n < nMax - t - u; n++){
                if (u == 0){
                    R.at(n)(t,u+1,0) = A(1)*R.at(n+1)(t,u,0);
                } else {
                    R.at(n)(t,u+1,0) = u*R.at(n+1)(t,u-1,0) + A(1)*R.at(n+1)(t,u,0);
                }
            }
        }
    }

    for(int t = 0; t < tMax+1; t++){
        for(int u = 0; u < uMax+1; u++){
            for(int v = 0; v < vMax; v++){
                for(int n = 0; n < nMax - t - u - v; n++){
                    if (v == 0){
                        R.at(n)(t,u,v+1) = A(2)*R.at(n+1)(t,u,v);
                    } else {
                        R.at(n)(t,u,v+1) = v*R.at(n+1)(t,u,v-1) + A(2)*R.at(n+1)(t,u,v);
                    }
                }
            }
        }
    }


//  Wrong!!
//    for (int t = 0; t < tMax; t++){
//        for (int n = 0; n < nMax - t; n++){
//            if (t == 0){
//                R.at(n)(t+1,0,0) = PC(0)*R.at(n+1)(t,0,0);
//            } else {
//                R.at(n)(t+1,0,0) = t*R.at(n+1)(t-1,0,0) + PC(0)*R.at(n+1)(t,0,0);
//            }
//        }
//    }

//    for (int u = 0; u < tMax; u++){
//        for (int n = 0; n < nMax - u; n++){
//            if (u == 0){
//                R.at(n)(0,u+1,0) = PC(1)*R.at(n+1)(0,u,0);
//            } else {
//                R.at(n)(0,u+1,0) = u*R.at(n+1)(0,u-1,0) + PC(1)*R.at(n+1)(0,u,0);
//            }
//        }
//    }

//    for (int u = 0; u < tMax+1; u++){
//        for (int v = 0; v < tMax; v++){
//            for (int n = 0; n < nMax - v; n++){
//                if (v == 0){
//                    R.at(n)(0,u,v+1) = PC(2)*R.at(n+1)(t,u,v);
//                } else {
//                    R.at(n)(0,u,v+1) = v*R.at(n+1)(t,u,v-1) + PC(2)*R.at(n+1)(t,u,v);
//                }
//            }
//        }
//    }



//    for(int tuv = 0; tuv < nMax+1; tuv++){
//        for(int n = 0; n < nMax; n++){
//            for(int t = 0; t < tMax; t++){
//                for(int u = 0; u < tMax; u++){
//                    for(int v = 0; v < tMax; v++){
//                        if(t + u + v == tuv){

//                            if(t == 0){
//                                R.at(n)(t+1,u,v) = PC(0)*R.at(n+1)(t,u,v);
//                            } else {
//                                R.at(n)(t+1,u,v) = t*R.at(n+1)(t-1,u,v) + PC(0)*R.at(n+1)(t,u,v);
//                            }

//                            if(u == 0){
//                                R.at(n)(t,u+1,v) = PC(1)*R.at(n+1)(t,u,v);
//                            } else {
//                                R.at(n)(t,u+1,v) = u*R.at(n+1)(t,u-1,v) + PC(1)*R.at(n+1)(t,u,v);
//                            }

//                            if(v == 0){
//                                R.at(n)(t,u,v+1) = PC(2)*R.at(n+1)(t,u,v);
//                            } else {
//                                R.at(n)(t,u,v+1) = v*R.at(n+1)(t,u,v-1) + PC(2)*R.at(n+1)(t,u,v);
//                            }

//                        }
//                    }
//                }
//            }
//        }
//    }

}

double Integrator::overlap(int i, int j, int k, int l, int m, int n)
{
    double p = alpha + beta;
    return E_AB[0](i,j,0)*E_AB[1](k,l,0)*E_AB[2](m,n,0)*pow(M_PI/p, 1.5);
}

double Integrator::kinetic(int i, int j, int k, int l, int m, int n)
{
    double Tij, Tkl, Tmn, kinetic;
    double p = alpha + beta;

    if (j < 2){
        Tij = 4*beta*beta*E_AB[0](i,j+2,0) - 2*beta*(2*j + 1)*E_AB[0](i,j,0);
    } else {
        Tij = 4*beta*beta*E_AB[0](i,j+2,0) - 2*beta*(2*j + 1)*E_AB[0](i,j,0) + j*(j-1)*E_AB[0](i,j-2,0);
    }

    if (l < 2){
        Tkl = 4*beta*beta*E_AB[1](k,l+2,0) - 2*beta*(2*l + 1)*E_AB[1](k,l,0);
    } else {
        Tkl = 4*beta*beta*E_AB[1](k,l+2,0) - 2*beta*(2*l + 1)*E_AB[1](k,l,0) + l*(l-1)*E_AB[1](k,l-2,0);
    }

    if (n < 2){
        Tmn = 4*beta*beta*E_AB[2](m,n+2,0) - 2*beta*(2*n + 1)*E_AB[2](m,n,0);
    } else {
        Tmn = 4*beta*beta*E_AB[2](m,n+2,0) - 2*beta*(2*n + 1)*E_AB[2](m,n,0) + n*(n-1)*E_AB[2](m,n-2,0);
    }

    kinetic = -0.5*(Tij*E_AB[1](k,l,0)*E_AB[2](m,n,0) + E_AB[0](i,j,0)*Tkl*E_AB[2](m,n,0) + E_AB[0](i,j,0)*E_AB[1](k,l,0)*Tmn)*pow(M_PI/p, 1.5);

    return kinetic;
}


double Integrator::coulomb1(int i, int j, int k, int l, int m, int n)
{
    int tMax = i + j;
    int uMax = k + l;
    int vMax = m + n;

    double p = alpha + beta;
    rowvec3 P = (alpha*Rnuclei.row(0) + beta*Rnuclei.row(1))/p;
    rowvec3 PC = P - Rnuclei.row(2); // tests coulomb1 for the particular case 1/rc = 1/Rnuclei.row(2)

    setR(p, PC, tMax, uMax, vMax);

    double value = 0;

    for (int t = 0; t < tMax + 1; t++){
        for (int u = 0; u < uMax + 1; u++){
            for (int v = 0; v < vMax + 1; v++){
                value += E_AB[0](i,j,t)*E_AB[1](k,l,u)*E_AB[2](m,n,v)*R.at(0)(t,u,v);
            }
        }
    }

    value *= 2*M_PI/p;

    return value;
}


double Integrator::coulomb2(int i1, int j1, int k1, int l1, int m1, int n1, int i2, int j2, int k2, int l2, int m2, int n2)
{
    int tMax = i1 + j1;
    int uMax = k1 + l1;
    int vMax = m1 + n1;
    int tauMax = i2 + j2;
    int nyMax = k2 + l2;
    int phiMax = m2 + n2;

    double p = alpha + beta;
    double q = gamma + delta;
    rowvec3 P = (alpha*Rnuclei.row(0) + beta*Rnuclei.row(1))/p;
    rowvec3 Q = (gamma*Rnuclei.row(2) + delta*Rnuclei.row(3))/q;
    rowvec PQ = P - Q;
    double a = p*q/(p + q);

    setR(a, PQ, tMax+tauMax, uMax+nyMax, vMax+phiMax);

    double value = 0;

    for (int t = 0; t < tMax+1; t++){
        for (int u = 0; u < uMax+1; u++){
            for (int v = 0; v < vMax+1; v++){
                for (int tau = 0; tau < tauMax+1; tau++){
                    for (int ny = 0; ny < nyMax+1; ny++){
                        for (int phi = 0; phi < phiMax+1; phi++){
                            value += (1 - 2*((tau + ny + phi)%2))*E_AB[0](i1,j1,t)*E_AB[1](k1,l1,u)*E_AB[2](m1,n1,v)
                                     *E_CD[0](i2,j2,tau)*E_CD[1](k2,l2,ny)*E_CD[2](m2,n2,phi)*R.at(0)(t+tau,u+ny,v+phi);
                        }
                    }
                }
            }
        }
    }

    value *= 2*pow(M_PI, 2.5)/(p*q*sqrt(p + q));

    return value;
}
