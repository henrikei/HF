#include <iostream>
#include <ctime>
#include <time.h>
#include <iomanip>
#include <fstream>
#include <libconfig.h++>
#include "hartreefock.h"
#include "rhf.h"
#include "uhf.h"
#include "system.h"
#include "integrator.h"
#include "boysfunction.h"
#include "basishandler.h"
#include "basisfunctions/basisfunctions.h"


using namespace std;
using namespace arma;
using namespace libconfig;

int main()
{
    Config cfg;
    cfg.readFile("../inFiles/configFiles/CH4_631Gs.cfg");
    string name = cfg.lookup("name");

    // Initialise basis functions and add them to basisHandler
    Setting& root = cfg.getRoot();
    Setting& atoms = root["atoms"];
    int nAtoms = atoms.getLength();
    rowvec position = zeros<rowvec>(3);
    mat nucleiPositions = zeros<mat>(nAtoms, 3);
    rowvec charges = zeros<rowvec>(nAtoms);
    string basis;
    BasisFunctions* basisfunctions;
    BasisHandler* basisHandler = new BasisHandler;
    for(int i = 0; i < nAtoms; i++){
        Setting& atom = atoms[i];
        atom.lookupValue("posX", position(0));
        atom.lookupValue("posY", position(1));
        atom.lookupValue("posZ", position(2));
        nucleiPositions.row(i) = position;
        atom.lookupValue("charge", charges(i));
        atom.lookupValue("basis", basis);
        basisfunctions = new BasisFunctions("../inFiles/basisSets/"+basis);
        basisfunctions->setPosition(position);
        basisHandler->addBasisFunctions(basisfunctions);
    }

    // Initialise system
    int nElectrons = cfg.lookup("nElectrons");
    System* system = new System(basisHandler, nucleiPositions, charges, nElectrons);

    // Solver type (RHF/UHF) and perturbation theory (true/false)
    HartreeFock* solver;
    string method = cfg.lookup("method");
    int pert = 1;
    if (cfg.lookup("perturbation")){
        pert = 2;}
    if (method == "RHF"){
        solver = new RHF(system, pert);
    } else if (method == "UHF"){
        solver = new UHF(system, pert);
    } else {
        cout << "Invalid choice of method" << endl;
        exit(EXIT_FAILURE);
    }
    solver->solve();
    cout << name << endl;
    cout << "Energy: " << setprecision(10) << solver->getEnergy() << endl;
    if (pert == 2){
        cout << "MP2 correlation energy: " << solver->getEnergyMP2() << endl;
    }



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

//    BasisFunctions* basis = new BasisFunctions("../inFiles/basisSets/O_431G.dat");
//    basis->setPosition(posO);
//    basisHandler->addBasisFunctions(basis);

//    basis = new BasisFunctions("../inFiles/basisSets/H_431G.dat");
//    basis->setPosition(posH1);
//    basisHandler->addBasisFunctions(basis);

//    basis = new BasisFunctions("../inFiles/basisSets/H_431G.dat");
//    basis->setPosition(posH2);
//    basisHandler->addBasisFunctions(basis);


//    System *system;
//    system = new System(basisHandler, nucleiPositions, charges, nElectrons);

//    RHF solver(system);
//    solver.solve();
//    cout << "Energy: " << solver.getEnergy() << endl;

//    clock_t end = clock();
//    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;



//    clock_t begin = clock();

//    double d = 1.4;
//    rowvec posA = {-d/2, 0.0, 0.0};
//    rowvec posB = {d/2, 0.0, 0.0};
//    rowvec charges = {1.0, 1.0};
//    int nElectrons = 2;

//    mat nucleiPositions = zeros<mat>(2,3);
//    nucleiPositions.row(0) = posA;
//    nucleiPositions.row(1) = posB;

//    BasisHandler* basisHandler = new BasisHandler;

//    BasisFunctions* basis = new BasisFunctions("../inFiles/basisSets/H_431G.dat");
//    basis->setPosition(posA);
//    basisHandler->addBasisFunctions(basis);

//    basis = new BasisFunctions("../inFiles/basisSets/H_431G.dat");
//    basis->setPosition(posB);
//    basisHandler->addBasisFunctions(basis);

//    System *system;
//    system = new System(basisHandler, nucleiPositions, charges, nElectrons);

//    RHF solver(system,2);
//    solver.solve();
//    cout << "Energy: " << setprecision(9) << solver.getEnergy() << endl;
//    cout << "MP2 correction: " << setprecision(9) << solver.getEnergyMP2() << endl;

//    clock_t end = clock();
//    cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;

    return 0;
}

