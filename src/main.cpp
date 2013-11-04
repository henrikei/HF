#include <iostream>
#include <iomanip>
#include <fstream>
#include "hartreefock.h"
#include "system.h"
#include "integrator.h"
#include "boysfunction.h"
#include "basishandler.h"
#include "basisfunctions/basisfunctions.h"
#include "basisfunctions/h_321g.h"
#include "basisfunctions/h_theijssen.h"

using namespace std;

int main()
{
    rowvec posA = {-0.5, 0.0, 0.0};
    rowvec posB = {0.5, 0.0, 0.0};

    mat nucleiPositions = zeros<mat>(2,3);
    nucleiPositions.row(0) = posA;
    nucleiPositions.row(1) = posB;

    BasisHandler* basisHandler = new BasisHandler;

    BasisFunctions* basis;
    basis = new H_Theijssen;
    basis->setPosition(posA);
    basisHandler->addBasisFunctions(basis);

    basis = new H_Theijssen;
    basis->setPosition(posB);
    basisHandler->addBasisFunctions(basis);


    System *system;
    system = new System(basisHandler, nucleiPositions);

    cout << basisHandler->getTotalNumOfBasisFunc() << endl;
    cout << basisHandler->getExponents(6) << endl;

    cout << system->getOneElectronIntegrals(0,0) << endl;
    cout << system->getTwoElectronIntegral(0,0,0,0) << endl;

    HartreeFock solver(system);
    solver.solve();
    cout << solver.getEnergy() << endl;


//    Integrator integrator;

//    rowvec3 RA = {1.2, 2.3, 3.4};
//    rowvec3 RB = {-1.3, 1.4, -2.4};
//    rowvec3 RC = {2.3, 0.9, 3.2};
//    rowvec3 RD = {-2.0, 1.9, 2.2};
//    double alpha = 0.2;
//    double beta = 0.3;
//    double gamma = 0.4;
//    double delta = 0.5;
//    int angMomMax = 3;

//    integrator.setPositionA(RA);
//    integrator.setPositionB(RB);
//    integrator.setPositionC(RC);
//    integrator.setPositionD(RD);
//    integrator.setAlpha(alpha);
//    integrator.setBeta(beta);
//    integrator.setGamma(gamma);
//    integrator.setDelta(delta);
//    integrator.setMaxAngMom(angMomMax);
//    integrator.setE();

//    cout << integrator.coulomb2(1,1,2,0,0,0,0,0,0,0,0,3) << endl;


//    mat nucleiPos = zeros<mat>(3,2);
//    nucleiPos(0,0) = -0.5;
//    nucleiPos(0,1) = 0.5;
//    HartreeFock H;
//    H.setPositions(nucleiPos);
//    H.solve();
//    cout << H.getOneElectronIntegral(0,0) << endl;
//    cout << H.getTwoElectronIntegral(0,0,0,0) << endl;
//    cout << H.getEnergy() << endl;




//    // Setup
//    HartreeFock H;
//    mat nucleiPos = zeros<mat>(3,2);
//    double distMin = 0.1;
//    double distMax = 6.0;
//    int nStep = 100;

//    ofstream out;

//    // Calculation

//    double distStep = (distMax - distMin)/(nStep-1);
//    out.open("../Results/hydrogen_energy_sto4g.dat");

//    for (int i = 0; i < nStep; i++){

//        nucleiPos(0,0) = -(distMin + i*distStep)/2;
//        nucleiPos(0,1) = (distMin + i*distStep)/2;

//        H.setPositions(nucleiPos);
//        H.solve();

//        out << distMin + i*distStep << "   " << H.getEnergy() << endl;
//    }
//    out.close();

    return 0;
}

