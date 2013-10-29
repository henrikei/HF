#include <iostream>
#include <iomanip>
#include <fstream>
#include "hartreefock.h"
#include "system.h"
#include "primitivebasis.h"
#include "integrator.h"
#include "boysfunction.h"


using namespace std;

int main()
{
    Integrator integrator;
    rowvec3 RA = {1.2, 2.3, 3.4};
    integrator.setPositionA(RA);
    rowvec3 RB = {-1.3, 1.4, -2.4};
    integrator.setPositionB(RB);
    rowvec3 RC = {2.3, 0.9, 3.2};
    integrator.setAlpha(0.2);
    integrator.setBeta(0.3);
    integrator.setMaxAngMom(2);
    integrator.setE();
    integrator.setR(RC);

    cout << integrator.R.at(0) << endl;
    cout << integrator.coulomb1(1,1,1,1,0,0) << endl;


//    mat nucleiPos = zeros<mat>(3,2);
//    nucleiPos(0,0) = -0.5;
//    nucleiPos(0,1) = 0.5;
//    HartreeFock H;
//    H.setPositions(nucleiPos);
//    H.solve();
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

