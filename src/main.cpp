#include <iostream>
#include <time.h>
#include <iomanip>
#include <fstream>
#include "hartreefock.h"
#include "system.h"
#include "integrator.h"
#include "boysfunction.h"
#include "basishandler.h"
#include "basisfunctions/basisfunctions.h"
#include "basisfunctions/h_321g.h"
#include "basisfunctions/h_431g.h"
#include "basisfunctions/h_theijssen.h"
#include "basisfunctions/o_321g.h"
#include "basisfunctions/o_431g.h"
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

//    HartreeFock solver(system);
//    solver.solve();
//    cout << "Energy: " << solver.getEnergy() << endl;

//    clock_t end = clock();
//    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;





//    clock_t begin = clock();

//    double d = 2.28;
//    rowvec posA = {-d/2, 0.0, 0.0};
//    rowvec posB = {d/2, 0.0, 0.0};
//    rowvec charges = {8.0, 8.0};
//    int nElectrons = 16;

//    mat nucleiPositions = zeros<mat>(2,3);
//    nucleiPositions.row(0) = posA;
//    nucleiPositions.row(1) = posB;

//    BasisHandler* basisHandler = new BasisHandler;

//    BasisFunctions* basis;
//    basis = new O_431G;
//    basis->setPosition(posA);
//    basisHandler->addBasisFunctions(basis);

//    basis = new O_431G;
//    basis->setPosition(posB);
//    basisHandler->addBasisFunctions(basis);


//    System *system;
//    system = new System(basisHandler, nucleiPositions, charges, nElectrons);

//    HartreeFock solver(system);
//    solver.solve();
//    cout << "Energy: " << solver.getEnergy() << endl;

//    clock_t end = clock();
//    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;


    clock_t begin = clock();

    rowvec posA = {0.0, 0.0, 0.0};
    rowvec charges = {8.0};
    int nElectrons = 8;

    mat nucleiPositions = zeros<mat>(1,3);
    nucleiPositions.row(0) = posA;

    BasisHandler* basisHandler = new BasisHandler;

    BasisFunctions* basis;
    basis = new O_431G;
    basis->setPosition(posA);
    basisHandler->addBasisFunctions(basis);

    System *system;
    system = new System(basisHandler, nucleiPositions, charges, nElectrons);

    HartreeFock solver(system);
    solver.solve();
    cout << "Energy: " << solver.getEnergy() << endl;

    clock_t end = clock();
    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;

    return 0;
}

