#include "integrator.h"

Integrator::Integrator(int angMomMax)
{
    m_boys = new BoysFunction(angMomMax);
    m_nucleusPosition = zeros<rowvec>(3);
    m_angMom = angMomMax;
}

Integrator::~Integrator()
{
    delete m_boys;
}

void Integrator::setNucleusPosition(rowvec3 RC)
{
    m_nucleusPosition = RC;
}

void Integrator::setPrimitiveA(Primitive *primitiveA)
{
    m_primitiveA = primitiveA;
}

void Integrator::setPrimitiveB(Primitive *primitiveB)
{
    m_primitiveB = primitiveB;
}

void Integrator::setPrimitiveC(Primitive *primitiveC)
{
    m_primitiveC = primitiveC;
}

void Integrator::setPrimitiveD(Primitive *primitiveD)
{
    m_primitiveD = primitiveD;
}

void Integrator::setE_AB(string intType)
{
    setE(m_E_AB, m_primitiveA, m_primitiveB, intType);
}



void Integrator::setE_CD(string intType)
{
    setE(m_E_CD, m_primitiveC, m_primitiveD, intType);
}



void Integrator::setR(double a, rowvec3 A, int tMax, int uMax, int vMax)
{
    int nMax = tMax + uMax + vMax;

    m_boys->setx(a*dot(A,A));

    cube Rinit = zeros<cube>(tMax+1, uMax+1, vMax+1);

    // Generate initial Rs
    m_R.clear();
    for(int n = 0; n < nMax+1; n++){
        Rinit(0,0,0) = pow(-2*a, n)*m_boys->returnValue(n);
        m_R.push_back(Rinit);
    }


    for(int t = 0; t < tMax; t++){
        for(int n = 0; n < nMax - t; n++){
            if (t == 0){
                m_R.at(n)(t+1,0,0) = A(0)*m_R.at(n+1)(t,0,0);
            } else {
                m_R.at(n)(t+1,0,0) = t*m_R.at(n+1)(t-1,0,0) + A(0)*m_R.at(n+1)(t,0,0);
            }
        }
    }

    for(int t = 0; t < tMax+1; t++){
        for(int u = 0; u < uMax; u++){
            for(int n = 0; n < nMax - t - u; n++){
                if (u == 0){
                    m_R.at(n)(t,u+1,0) = A(1)*m_R.at(n+1)(t,u,0);
                } else {
                    m_R.at(n)(t,u+1,0) = u*m_R.at(n+1)(t,u-1,0) + A(1)*m_R.at(n+1)(t,u,0);
                }
            }
        }
    }

    for(int t = 0; t < tMax+1; t++){
        for(int u = 0; u < uMax+1; u++){
            for(int v = 0; v < vMax; v++){
                for(int n = 0; n < nMax - t - u - v; n++){
                    if (v == 0){
                        m_R.at(n)(t,u,v+1) = A(2)*m_R.at(n+1)(t,u,v);
                    } else {
                        m_R.at(n)(t,u,v+1) = v*m_R.at(n+1)(t,u,v-1) + A(2)*m_R.at(n+1)(t,u,v);
                    }
                }
            }
        }
    }
}

double Integrator::overlap()
{
    int i = m_primitiveA->getPow()(0);
    int j = m_primitiveB->getPow()(0);
    int k = m_primitiveA->getPow()(1);
    int l = m_primitiveB->getPow()(1);
    int m = m_primitiveA->getPow()(2);
    int n = m_primitiveB->getPow()(2);

    return overlap(i, j, k, l, m, n); // Parameters other than i,j,k... are set inside overlap(i,j,k,l,m,n)
}

double Integrator::kinetic()
{
    int i = m_primitiveA->getPow()(0);
    int j = m_primitiveB->getPow()(0);
    int k = m_primitiveA->getPow()(1);
    int l = m_primitiveB->getPow()(1);
    int m = m_primitiveA->getPow()(2);
    int n = m_primitiveB->getPow()(2);

    return kinetic(i, j, k, l, m, n); // Parameters other than i,j,k... are set inside kinetic(i,j,k,l,m,n)
}

double Integrator::coulomb1()
{
    int i = m_primitiveA->getPow()(0);
    int j = m_primitiveB->getPow()(0);
    int k = m_primitiveA->getPow()(1);
    int l = m_primitiveB->getPow()(1);
    int m = m_primitiveA->getPow()(2);
    int n = m_primitiveB->getPow()(2);

    return coulomb1(i, j, k, l, m, n); // Parameters other than i,j,k... are set inside kinetic(i,j,k,l,m,n)
}

double Integrator::coulomb2()
{
    int i1 = m_primitiveA->getPow()(0);
    int j1 = m_primitiveB->getPow()(0);
    int k1 = m_primitiveA->getPow()(1);
    int l1 = m_primitiveB->getPow()(1);
    int m1 = m_primitiveA->getPow()(2);
    int n1 = m_primitiveB->getPow()(2);
    int i2 = m_primitiveC->getPow()(0);
    int j2 = m_primitiveD->getPow()(0);
    int k2 = m_primitiveC->getPow()(1);
    int l2 = m_primitiveD->getPow()(1);
    int m2 = m_primitiveC->getPow()(2);
    int n2 = m_primitiveD->getPow()(2);

    int tMax = i1 + j1;
    int uMax = k1 + l1;
    int vMax = m1 + n1;
    int tauMax = i2 + j2;
    int nyMax = k2 + l2;
    int phiMax = m2 + n2;

    double alpha = m_primitiveA->getExp();
    double beta = m_primitiveB->getExp();
    double gamma = m_primitiveC->getExp();
    double delta = m_primitiveD->getExp();

    rowvec A = m_primitiveA->getPos();
    rowvec B = m_primitiveB->getPos();
    rowvec C = m_primitiveC->getPos();
    rowvec D = m_primitiveD->getPos();

    double p = alpha + beta;
    double q = gamma + delta;
    rowvec3 P = (alpha*A + beta*B)/p;
    rowvec3 Q = (gamma*C + delta*D)/q;
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
                            value += (1 - 2*((tau + ny + phi)%2))*m_E_AB[0](i1,j1,t)*m_E_AB[1](k1,l1,u)*m_E_AB[2](m1,n1,v)
                                    *m_E_CD[0](i2,j2,tau)*m_E_CD[1](k2,l2,ny)*m_E_CD[2](m2,n2,phi)*m_R.at(0)(t+tau,u+ny,v+phi);
                        }
                    }
                }
            }
        }
    }

    value *= 2*pow(M_PI, 2.5)/(p*q*sqrt(p + q));

    return value;
}

double Integrator::overlap(int i, int j, int k, int l, int m, int n)
{
    double p = m_primitiveA->getExp() + m_primitiveB->getExp();
    return m_E_AB[0](i,j,0)*m_E_AB[1](k,l,0)*m_E_AB[2](m,n,0)*pow(M_PI/p, 1.5);
}

double Integrator::kinetic(int i, int j, int k, int l, int m, int n)
{
    double Tij, Tkl, Tmn, kinetic;
    double alpha = m_primitiveA->getExp();
    double beta = m_primitiveB->getExp();
    double p = alpha + beta;

    if (j < 2){
        Tij = 4*beta*beta*m_E_AB[0](i,j+2,0) - 2*beta*(2*j + 1)*m_E_AB[0](i,j,0);
    } else {
        Tij = 4*beta*beta*m_E_AB[0](i,j+2,0) - 2*beta*(2*j + 1)*m_E_AB[0](i,j,0) + j*(j-1)*m_E_AB[0](i,j-2,0);
    }

    if (l < 2){
        Tkl = 4*beta*beta*m_E_AB[1](k,l+2,0) - 2*beta*(2*l + 1)*m_E_AB[1](k,l,0);
    } else {
        Tkl = 4*beta*beta*m_E_AB[1](k,l+2,0) - 2*beta*(2*l + 1)*m_E_AB[1](k,l,0) + l*(l-1)*m_E_AB[1](k,l-2,0);
    }

    if (n < 2){
        Tmn = 4*beta*beta*m_E_AB[2](m,n+2,0) - 2*beta*(2*n + 1)*m_E_AB[2](m,n,0);
    } else {
        Tmn = 4*beta*beta*m_E_AB[2](m,n+2,0) - 2*beta*(2*n + 1)*m_E_AB[2](m,n,0) + n*(n-1)*m_E_AB[2](m,n-2,0);
    }

    kinetic = -0.5*(Tij*m_E_AB[1](k,l,0)*m_E_AB[2](m,n,0) + m_E_AB[0](i,j,0)*Tkl*m_E_AB[2](m,n,0) + m_E_AB[0](i,j,0)*m_E_AB[1](k,l,0)*Tmn)*pow(M_PI/p, 1.5);

    return kinetic;
}

double Integrator::coulomb1(int i, int j, int k, int l, int m, int n)
{
    int tMax = i + j;
    int uMax = k + l;
    int vMax = m + n;

    rowvec3 A = m_primitiveA->getPos();
    rowvec3 B = m_primitiveB->getPos();

    double alpha = m_primitiveA->getExp();
    double beta = m_primitiveB->getExp();
    double p = alpha + beta;
    rowvec3 P = (alpha*A + beta*B)/p;
    rowvec3 PC = P - m_nucleusPosition; // tests coulomb1 for the particular case 1/rc = 1/Rnuclei.row(2)

    setR(p, PC, tMax, uMax, vMax);

    double value = 0;

    for (int t = 0; t < tMax + 1; t++){
        for (int u = 0; u < uMax + 1; u++){
            for (int v = 0; v < vMax + 1; v++){
                value += m_E_AB[0](i,j,t)*m_E_AB[1](k,l,u)*m_E_AB[2](m,n,v)*m_R.at(0)(t,u,v);
            }
        }
    }

    value *= 2*M_PI/p;

    return value;
}

void Integrator::setE(cube E[], Primitive *primitive1, Primitive *primitive2, string intType)
{
    double alpha = primitive1->getExp();
    double beta = primitive2->getExp();

    rowvec3 A = primitive1->getPos();
    rowvec3 B = primitive2->getPos();

    double p = alpha + beta;
    rowvec3 P = (alpha*A + beta*B)/p;
    rowvec3 AB = A - B;
    rowvec3 PA = P - A;
    rowvec3 PB = P - B;

    E[0] = zeros<cube>(m_angMom+3,m_angMom+3,2*(m_angMom+2)+2);
    E[1] = zeros<cube>(m_angMom+3,m_angMom+3,2*(m_angMom+2)+2);
    E[2] = zeros<cube>(m_angMom+3,m_angMom+3,2*(m_angMom+2)+2);

    // Loop over x-, y- and z- coordinates

    for(int dir = 0; dir < 3; dir++){

        int iMax;
        if (intType == "oneParticle"){
            iMax = m_angMom + 2;
        } else if (intType == "twoParticle"){
            iMax = m_angMom;
        } else {
            cout << "Error: Integration type for E-coeffs not specified" << endl;
            exit(EXIT_FAILURE);
        }

        // First loop over t and i with j = 0

        E[dir](0,0,0) = exp(-alpha*beta*AB(dir)*AB(dir)/p);
        for (int i = 0; i < iMax; i++){
            // Treat the case t=0 separately due to (t-1) term
            E[dir](i+1,0,0) = PA(dir)*E[dir](i,0,0) + E[dir](i,0,1);
            for (int t = 1; t <= i + 0 + 1; t++){
                E[dir](i+1,0,t) = E[dir](i,0,t-1)/(2*p) + PA(dir)*E[dir](i,0,t) + (t+1)*E[dir](i,0,t+1);
            }
        }

        // Second loop over t and j and i

        // Must here let i <= l because the forward loop is on index j
        for (int i = 0; i <= iMax; i++){
            for (int j = 0; j < iMax; j++){
                E[dir](i,j+1,0) = PB(dir)*E[dir](i,j,0) + E[dir](i,j,1);
                for (int t = 1; t <= i + j + 1; t++){
                    E[dir](i,j+1,t) = E[dir](i,j,t-1)/(2*p) + PB(dir)*E[dir](i,j,t) + (t+1)*E[dir](i,j,t+1);
                }
            }
        }
    }
}
