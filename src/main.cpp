#include <iostream>
#include <ctime>
#include <time.h>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <libconfig.h++>
#include "hartreefock.h"
#include "rhf.h"
#include "uhf.h"
#include "system.h"
#include "integrator.h"
#include "boysfunction.h"
#include "basishandler.h"
#include "basisfunctions/basisfunctions.h"
#include "basisfunctions/basisfunctions2.h"
#include "basisfunctions/contracted.h"
#include "basisfunctions/primitive.h"


using namespace std;
using namespace arma;
using namespace libconfig;

int main()
{
    //    BasisFunctions2 basis;
    //    rowvec3 position = {1, 2, 3};
    //    basis.addContracteds("../inFiles/basisSets/O_431G.dat", position);

    //    int M = basis.getNumOfContracteds();
    //    for (int i = 0; i < M; i++){
    //        Contracted *contracted = basis.getContracted(i);
    //        int L = contracted->getNumOfPrimitives();
    //        cout << "Contracted " << i+1 << " :" << endl;
    //        for (int j = 0; j < L; j++){
    //            Primitive *primitive = contracted->getPrimitive(j);
    //            cout << primitive->getExp() << " , " << primitive->getCoeff() << " , " << primitive->getPos() << " , " << primitive->getPow() << endl << endl;
    //        }
    //    }

//    clock_t begin = clock();

//    double d = 1.1795265999544056;
//    rowvec posC = {0.0, 0.0, 0.0};
//    rowvec posH1 = {d, d, d};
//    rowvec posH2 = {-d, -d, d};
//    rowvec posH3 = {d, -d, -d};
//    rowvec posH4 = {-d, d, -d};
//    rowvec charges = {6.0, 1.0, 1.0, 1.0, 1.0};
//    int nElectrons = 10;

//    mat nucleiPositions = zeros<mat>(5,3);
//    nucleiPositions.row(0) = posC;
//    nucleiPositions.row(1) = posH1;
//    nucleiPositions.row(2) = posH2;
//    nucleiPositions.row(3) = posH3;
//    nucleiPositions.row(4) = posH4;

//    BasisFunctions2 *basisFunctions = new BasisFunctions2;
//    basisFunctions->addContracteds("../inFiles/basisSets/C_631Gs.dat", posC);
//    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", posH1);
//    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", posH2);
//    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", posH3);
//    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", posH4);

//    System *system;
//    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

//    RHF solver(system);
//    solver.solve();
//    cout << "Energy: " << setprecision(10) <<  solver.getEnergy() << endl;

//    clock_t end = clock();
//    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;



//    clock_t begin = clock();

//    double d = 2.028;//2.0262;
//    double s = sin(2*M_PI/3);
//    double c = cos(2*M_PI/3);
//    rowvec posC = {0.0, 0.0, 0.0};
//    rowvec posH1 = {d, 0.0, 0.0};
//    rowvec posH2 = {c*d, s*d, 0.0};
//    rowvec posH3 = {c*d, -s*d, 0.0};
//    rowvec charges = {6.0, 1.0, 1.0, 1.0};
//    int nElectrons = 9;

//    mat nucleiPositions = zeros<mat>(4,3);
//    nucleiPositions.row(0) = posC;
//    nucleiPositions.row(1) = posH1;
//    nucleiPositions.row(2) = posH2;
//    nucleiPositions.row(3) = posH3;

//    BasisFunctions2 *basisFunctions = new BasisFunctions2;
//    basisFunctions->addContracteds("../inFiles/basisSets/C_631Gs.dat", posC);
//    basisFunctions->addContracteds("../inFiles/basisSets/H_631G.dat", posH1);
//    basisFunctions->addContracteds("../inFiles/basisSets/H_631G.dat", posH2);
//    basisFunctions->addContracteds("../inFiles/basisSets/H_631G.dat", posH3);

//    System *system;
//    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

//    UHF solver(system,2);
//    solver.solve();
//    cout << "Energy: " << setprecision(7) <<  solver.getEnergy() << endl;

//    clock_t end = clock();
//    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;




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

//    BasisFunctions2* basisFunctions = new BasisFunctions2;
//    basisFunctions->addContracteds("../inFiles/basisSets/O_431G.dat", posO);
//    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", posH1);
//    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", posH2);

//    System *system;
//    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

//    RHF solver(system);
//    solver.solve();
//    cout << "Energy: " << solver.getEnergy() << endl;

//    clock_t end = clock();
//    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;




    clock_t begin = clock();

    double d = 1.889725989;
    rowvec posH = {-0.5*d, 0.0, 0.0};
    rowvec posF = {0.5*d, 0.0, 0.0};
    rowvec charges = {1.0, 9.0};
    int nElectrons = 10;

    mat nucleiPositions = zeros<mat>(2,3);
    nucleiPositions.row(0) = posH;
    nucleiPositions.row(1) = posF;

    BasisFunctions2* basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", posH);
    basisFunctions->addContracteds("../inFiles/basisSets/F_631Gs.dat", posF);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    RHF solver(system,2);
    solver.solve();
    cout << "Energy: " << setprecision(9) << solver.getEnergy() << endl;

    clock_t end = clock();
    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;

    return 0;
}

