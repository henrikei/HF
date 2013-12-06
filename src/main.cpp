#include <iostream>
#include <time.h>
#include <iomanip>
#include <fstream>
#include "hartreefock.h"
#include "rhf.h"
#include "uhf.h"
#include "system.h"
#include "integrator.h"
#include "boysfunction.h"
#include "basishandler.h"
#include "basisfunctions/basisfunctions.h"
#include "basisfunctions/h_321g.h"
#include "basisfunctions/h_431g.h"
#include "basisfunctions/h_theijssen.h"
#include "basisfunctions/h_6311gss.h"
#include "basisfunctions/o_321g.h"
#include "basisfunctions/o_431g.h"
#include "basisfunctions/o_6311g.h"
#include "basisfunctions/o_6311gss.h"
#include "basisfunctions/n_431g.h"

using namespace std;

int main()
{
//    clock_t begin = clock();

//    rowvec posO = {0.0, 0.0, 0.0};
//    rowvec posH1 = {1.797, 0.0, 0.0};
//    rowvec posH2 = {-0.448, 1.740, 0.0};
//    rowvec charges = {8.0, 1.0, 1.0};
//    int nElectrons = 10;

//    mat nucleiPositions = zeros<mat>(3,3);
//    nucleiPositions.row(0) = posO;
//    nucleiPositions.row(1) = posH1;
//    nucleiPositions.row(2) = posH2;

//    BasisHandler* basisHandler = new BasisHandler;

//    BasisFunctions* basis;
//    basis = new O_431G;
//    basis->setPosition(posO);
//    basisHandler->addBasisFunctions(basis);

//    basis = new H_431G;
//    basis->setPosition(posH1);
//    basisHandler->addBasisFunctions(basis);

//    basis = new H_431G;
//    basis->setPosition(posH2);
//    basisHandler->addBasisFunctions(basis);


//    System *system;
//    system = new System(basisHandler, nucleiPositions, charges, nElectrons);

//    RHF solver(system);
//    solver.solve();
//    cout << "Energy: " << solver.getEnergy() << endl;

//    clock_t end = clock();
//    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;




//    double dMin = 0.5;
//    double dMax = 6.0;
//    int nPoints = 200;
//    double dDelta = (dMax - dMin)/(double (nPoints-1));

//    fstream ofile;
//    ofile.open("../out/H2_UHF_6311Gss_PT.dat");
//    clock_t begin = clock();

//    for (int i = 0; i < nPoints; i++){
//        double d = dMin + i*dDelta;
//        rowvec posA = {-d/2, 0.0, 0.0};
//        rowvec posB = {d/2, 0.0, 0.0};
//        rowvec charges = {1.0, 1.0};
//        int nElectrons = 2;

//        mat nucleiPositions = zeros<mat>(2,3);
//        nucleiPositions.row(0) = posA;
//        nucleiPositions.row(1) = posB;

//        BasisHandler* basisHandler = new BasisHandler;

//        BasisFunctions* basis;
//        basis = new H_6311Gss;
//        basis->setPosition(posA);
//        basisHandler->addBasisFunctions(basis);

//        basis = new H_6311Gss;
//        basis->setPosition(posB);
//        basisHandler->addBasisFunctions(basis);


//        System *system;
//        system = new System(basisHandler, nucleiPositions, charges, nElectrons);

//        UHF solver(system);
//        solver.solve();
//        ofile << d << "  " << solver.getEnergy() << endl;

//        delete basisHandler;
//        delete basis;
//        delete system;
//    }
//    ofile.close();
//    clock_t end = clock();
//    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;




    clock_t begin = clock();

    double d = 2.050;
    rowvec posA = {-d/2, 0.0, 0.0};
    rowvec posB = {d/2, 0.0, 0.0};
    rowvec charges = {7.0, 7.0};
    int nElectrons = 14;

    mat nucleiPositions = zeros<mat>(2,3);
    nucleiPositions.row(0) = posA;
    nucleiPositions.row(1) = posB;

    BasisHandler* basisHandler = new BasisHandler;

    BasisFunctions* basis;
    basis = new N_431G;
    basis->setPosition(posA);
    basisHandler->addBasisFunctions(basis);

    basis = new N_431G;
    basis->setPosition(posB);
    basisHandler->addBasisFunctions(basis);

    System *system;
    system = new System(basisHandler, nucleiPositions, charges, nElectrons);

    UHF solver(system);
    solver.solve();
    cout << "Energy: " << solver.getEnergy() << endl;

    clock_t end = clock();
    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;

    return 0;
}

