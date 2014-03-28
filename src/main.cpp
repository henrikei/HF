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
#include "basisfunctions/basisfunctions2.h"
#include "basisfunctions/contracted.h"
#include "basisfunctions/primitive.h"
#include "minimizer/func.h"
#include "minimizer/minimizer.h"
#include "minimizer/twodimtest.h"
#include "minimizer/hartreefockfunc.h"


using namespace std;
using namespace arma;
using namespace libconfig;

int main()
{
    string run = "FCl";

    if (run == "CH4"){

        clock_t begin = clock();

        double d = 1.1795265999544056;
        rowvec posC = {0.0, 0.0, 0.0};
        rowvec posH1 = {d, d, d};
        rowvec posH2 = {-d, -d, d};
        rowvec posH3 = {d, -d, -d};
        rowvec posH4 = {-d, d, -d};
        rowvec charges = {6.0, 1.0, 1.0, 1.0, 1.0};
        int nElectrons = 10;

        mat nucleiPositions = zeros<mat>(5,3);
        nucleiPositions.row(0) = posC;
        nucleiPositions.row(1) = posH1;
        nucleiPositions.row(2) = posH2;
        nucleiPositions.row(3) = posH3;
        nucleiPositions.row(4) = posH4;

        BasisFunctions2 *basisFunctions = new BasisFunctions2;
        basisFunctions->addContracteds("../inFiles/basisSets/C_631Gs.dat", 0);
        basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 1);
        basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 2);
        basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 3);
        basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 4);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        RHF solver(system);
        solver.solve();
        cout << "Energy: " << setprecision(10) <<  solver.getEnergy() << endl;

        clock_t end = clock();
        cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;

    } else if (run == "CH3"){

        clock_t begin = clock();

        double d = 2.028;//2.0262;
        double s = sin(2*M_PI/3);
        double c = cos(2*M_PI/3);
        rowvec posC = {0.0, 0.0, 0.0};
        rowvec posH1 = {d, 0.0, 0.0};
        rowvec posH2 = {c*d, s*d, 0.0};
        rowvec posH3 = {c*d, -s*d, 0.0};
        rowvec charges = {6.0, 1.0, 1.0, 1.0};
        int nElectrons = 9;

        mat nucleiPositions = zeros<mat>(4,3);
        nucleiPositions.row(0) = posC;
        nucleiPositions.row(1) = posH1;
        nucleiPositions.row(2) = posH2;
        nucleiPositions.row(3) = posH3;

        BasisFunctions2 *basisFunctions = new BasisFunctions2;
        basisFunctions->addContracteds("../inFiles/basisSets/C_631Gs.dat", 0);
        basisFunctions->addContracteds("../inFiles/basisSets/H_631G.dat", 1);
        basisFunctions->addContracteds("../inFiles/basisSets/H_631G.dat", 2);
        basisFunctions->addContracteds("../inFiles/basisSets/H_631G.dat", 3);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        UHF solver(system,2);
        solver.solve();
        cout << "Energy: " << setprecision(7) <<  solver.getEnergy() << endl;

        clock_t end = clock();
        cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;

    } else if (run == "H2O"){

        clock_t begin = clock();

        rowvec posO = {0.0, 0.0, 0.0};
        rowvec posH1 = {1.797, 0.0, 0.0};
        rowvec posH2 = {-0.448, 1.740, 0.0};
        rowvec charges = {8.0, 1.0, 1.0};
        int nElectrons = 10;

        mat nucleiPositions = zeros<mat>(3,3);
        nucleiPositions.row(0) = posO;
        nucleiPositions.row(1) = posH1;
        nucleiPositions.row(2) = posH2;

        BasisFunctions2* basisFunctions = new BasisFunctions2;
        basisFunctions->addContracteds("../inFiles/basisSets/O_431G.dat", 0);
        basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 1);
        basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 2);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        RHF solver(system);
        solver.solve();
        cout << "Energy: " << solver.getEnergy() << endl;

        clock_t end = clock();
        cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;

    } else if (run =="H2"){

        clock_t begin = clock();

        double d = 1.4;
        rowvec posH1 = {-d/2, 0.0, 0.0};
        rowvec posH2 = {d/2, 0.0, 0.0};
        rowvec charges = {1.0, 1.0};
        int nElectrons = 2;

        mat nucleiPositions = zeros<mat>(2,3);
        nucleiPositions.row(0) = posH1;
        nucleiPositions.row(1) = posH2;

        BasisFunctions2* basisFunctions = new BasisFunctions2;
        basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 0);
        basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 1);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        UHF solver(system,3);
        solver.solve();
        cout << "Energy: " << solver.getEnergy() << endl;
        cout << "energyMP2: " << solver.getEnergyMP2() << endl;
        cout << "energyMP3: " << solver.getEnergyMP2()+solver.getEnergyMP3() << endl;

        clock_t end = clock();
        cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;

    } else if (run == "HF") {

        clock_t begin = clock();

        double d = 1.735524348;
        rowvec posH = {-0.5*d, 0.0, 0.0};
        rowvec posF= {0.5*d, 0.0, 0.0};
        rowvec charges = {1.0, 9.0};
        int nElectrons = 10;

        mat nucleiPositions = zeros<mat>(2,3);
        nucleiPositions.row(0) = posH;
        nucleiPositions.row(1) = posF;

        BasisFunctions2* basisFunctions = new BasisFunctions2;
        basisFunctions->addContracteds("../inFiles/basisSets/H_6311Gss.dat", 0);
        basisFunctions->addContracteds("../inFiles/basisSets/F_6311Gs.dat", 1);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        RHF solver(system,2);
        solver.solve();
        cout << "Energy: " << setprecision(9) << solver.getEnergy() << endl;

        clock_t end = clock();
        cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;

        delete system;
        delete basisFunctions;

    } else if (run == "FCl") {

        clock_t begin = clock();

        double d = 3.154519593;
        rowvec posF = {-0.5*d, 0.0, 0.0};
        rowvec posCl= {0.5*d, 0.0, 0.0};
        rowvec charges = {9.0, 17.0};
        int nElectrons = 26;

        mat nucleiPositions = zeros<mat>(2,3);
        nucleiPositions.row(0) = posF;
        nucleiPositions.row(1) = posCl;

        BasisFunctions2* basisFunctions = new BasisFunctions2;
        basisFunctions->addContracteds("../inFiles/basisSets/F_6311Gs.dat", 0);
        basisFunctions->addContracteds("../inFiles/basisSets/Cl_6311Gs.dat", 1);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        UHF solver(system,3);
        solver.solve();
        cout << "Energy: " << setprecision(9) << solver.getEnergy() << endl;

        clock_t end = clock();
        cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;

    } else if (run == "H2O_Minimize"){

        rowvec O = {0.0,0.0,0.0};
        rowvec H1 = {1.0,0.0,0.0};
        rowvec H2 = {0.0,1.0,0.0};
        rowvec charges = {8.0,1.0,1.0};
        int nElectrons = 10;
        mat nucleiPositions = zeros<mat>(3,3);
        nucleiPositions.row(0) = O;
        nucleiPositions.row(1) = H1;
        nucleiPositions.row(2) = H2;
        BasisFunctions2* basisFunctions = new BasisFunctions2;
        basisFunctions->addContracteds("../inFiles/basisSets/O_431G.dat", 0);
        basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 1);
        basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 2);
        System *system = new System(basisFunctions, nucleiPositions, charges, nElectrons);
        RHF *solver = new RHF(system);
        HartreeFockFunc *func = new HartreeFockFunc(solver, system);
        Minimizer *minimizer = new Minimizer(func);
        minimizer->solve();
        cout << system->getNucleiPositions() << endl;
        cout << minimizer->getMinValue() << endl;

    } else {
        cout << "No valid run selected." << endl;
    }


    return 0;
}

