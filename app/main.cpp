#include <iostream>
#include <ctime>
#include <time.h>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <libconfig.h++>
#ifdef RUN_MPI
#include <mpi.h>
#endif
#include <hartreefock/hartreefock.h>
#include <hartreefock/rhf.h>
#include <hartreefock/uhf.h>
#include <perturbation/rmp.h>
#include <perturbation/ump.h>
#include <system/system.h>
#include <integrator/integrator.h>
#include <boysfunction/boysfunction.h>
#include <basisfunctions/basisfunctions.h>
#include <basisfunctions/contracted.h>
#include <basisfunctions/primitive.h>
#include <minimizer/func.h>
#include <minimizer/minimizer.h>
#include <minimizer/twodimtest.h>
#include <minimizer/hartreefockfunc.h>
#include <density/density.h>


using namespace std;
using namespace arma;
using namespace libconfig;

//#define RUN_MPI

int main()
{
    // MPI
    int my_rank = 0;
#ifdef RUN_MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

    string run = "H2_potential";


    clock_t begin = clock();

    if (run == "CH4"){
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

        BasisFunctions *basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/C_631Gs.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631Gss.dat", 1);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631Gss.dat", 2);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631Gss.dat", 3);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631Gss.dat", 4);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        RHF solver(system);
        solver.solve();

        if (my_rank == 0){
            cout << "Energy: " << setprecision(10) <<  solver.getEnergy() << endl;
        }

    } else if (run == "NH4"){

        double d = 1.1150365517919136;
        rowvec posN = {0.0, 0.0, 0.0};
        rowvec posH1 = {d, d, d};
        rowvec posH2 = {-d, -d, d};
        rowvec posH3 = {d, -d, -d};
        rowvec posH4 = {-d, d, -d};
        rowvec charges = {7.0, 1.0, 1.0, 1.0, 1.0};
        int nElectrons = 10;

        mat nucleiPositions = zeros<mat>(5,3);
        nucleiPositions.row(0) = posN;
        nucleiPositions.row(1) = posH1;
        nucleiPositions.row(2) = posH2;
        nucleiPositions.row(3) = posH3;
        nucleiPositions.row(4) = posH4;

        BasisFunctions *basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/N_631++Gs.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631++Gss.dat", 1);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631++Gss.dat", 2);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631++Gss.dat", 3);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631++Gss.dat", 4);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        RMP solver(system,3,2);
        solver.solve();

        if (my_rank == 0){
            cout << "Energy: " << setprecision(10) <<  solver.getEnergyHF() + solver.getEnergy2order() + solver.getEnergy3order()<< endl;
        }

    } else if (run == "CH3"){

        double d = 2.04;//2.0262;
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

        BasisFunctions *basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/C_631Gs.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631G.dat", 1);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631G.dat", 2);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631G.dat", 3);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        UMP solver(system,2,2);
        solver.solve();

        if (my_rank == 0){
            cout << "Energy: " << setprecision(7) <<  solver.getEnergyHF() + solver.getEnergy2order() + solver.getEnergy3order() << endl;
        }

    } else if (run == "H2O"){

        rowvec posO = {0.0, 0.0, 0.0};
        rowvec posH1 = {1.797, 0.0, 0.0};
        rowvec posH2 = {-0.448, 1.740, 0.0};
        rowvec charges = {8.0, 1.0, 1.0};
        int nElectrons = 10;

        mat nucleiPositions = zeros<mat>(3,3);
        nucleiPositions.row(0) = posO;
        nucleiPositions.row(1) = posH1;
        nucleiPositions.row(2) = posH2;

        BasisFunctions* basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/O_6311++Gss.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_6311++Gss.dat", 1);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_6311++Gss.dat", 2);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        RHF solver(system);
        solver.solve();

        if (my_rank == 0){
            cout << "Energy: " << solver.getEnergy() << endl;
        }

    } else if (run =="H2"){

        double d = 2.8;
        rowvec posH1 = {-d/2, 0.0, 0.0};
        rowvec posH2 = {d/2, 0.0, 0.0};
        rowvec charges = {1.0, 1.0};
        int nElectrons = 2;

        mat nucleiPositions = zeros<mat>(2,3);
        nucleiPositions.row(0) = posH1;
        nucleiPositions.row(1) = posH2;

        BasisFunctions* basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631Gss.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631Gss.dat", 1);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        UHF solver(system);
        solver.solve();

        if (my_rank == 0){
            cout << "Energy: " << solver.getEnergy() << endl;

            field<mat> C = solver.getCoeff();
            rowvec R1 = {-6.0, -3.0, -3.0};
            rowvec R2 = {6.0, 3.0, 3.0};
            double dx = 8.0/100;
            Density density(basisFunctions, R1, R2, dx, dx, dx);
            density.printMolecularOrbital(C, 1,"../../Results/H2_potential/orbitals_d2p0/orbital2_UHF.dat");
        }

    } else if (run == "H2_potential") {

        rowvec charges = {1.0, 1.0};
        int nElectrons = 2;
        mat nucleiPositions = zeros<mat>(2,3);

        BasisFunctions* basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631Gss.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631Gss.dat", 1);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);
        UMP solver(system, 3);

        double dmin = 0.5;
        double dmax = 6.0;
        int nPoints = 551;
        double delta = (dmax - dmin)/(nPoints-1);
        double d, energy;
        rowvec3 posH1, posH2;

        ofstream ofile;
        ofile.open("../../Results/H2_potential/UHF_631Gss.dat");
        ofile << "Basis set: 6-31G**" << endl;
        ofile << "Distance   UHF   UMP2   UMP3" << endl;

        for (int i = 0; i < nPoints; i++){
            d = dmin + i*delta;
            posH1 = {-d/2, 0, 0}; posH2 = {d/2, 0, 0};
            nucleiPositions.row(0) = posH1; nucleiPositions.row(1) = posH2;
            system->setNucleiPositions(nucleiPositions);
            solver.solve();
            ofile << d << "  ";
            energy = solver.getEnergyHF();
            ofile << setprecision(10) << energy << "  ";
            energy += solver.getEnergy2order();
            ofile << setprecision(10) << energy << "  ";
            energy += solver.getEnergy3order();
            ofile << setprecision(10) << energy << endl;
        }

    } else if (run == "HF") {

        double d = 1.889725989*2;
        rowvec posH = {-0.5*d, 0.0, 0.0};
        rowvec posF= {0.5*d, 0.0, 0.0};
        rowvec charges = {1.0, 9.0};
        int nElectrons = 10;

        mat nucleiPositions = zeros<mat>(2,3);
        nucleiPositions.row(0) = posH;
        nucleiPositions.row(1) = posF;

        BasisFunctions* basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_631Gss.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/F_631Gs.dat", 1);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        UMP solver(system,2,2);
        solver.solve();

        if (my_rank == 0){
            cout << "Energy: " << setprecision(9) << solver.getEnergyHF() + solver.getEnergy2order() << endl;
        }

        delete system;
        delete basisFunctions;

    } else if (run =="O2") {

        double d = 2.281;//2.314158446;
        rowvec posO1 = {-0.5*d, 0.0, 0.0};
        rowvec posO2= {0.5*d, 0.0, 0.0};
        rowvec charges = {8.0, 8.0};
        int nElectrons = 16;

        mat nucleiPositions = zeros<mat>(2,3);
        nucleiPositions.row(0) = posO1;
        nucleiPositions.row(1) = posO2;

        BasisFunctions* basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/O_631Gs.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/O_631Gs.dat", 1);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        UHF solver(system);
        solver.solve();

        if (my_rank == 0){
            cout << "Energy: " << setprecision(9) << solver.getEnergy() << endl;
        }

        delete system;
        delete basisFunctions;

    } else if (run == "N2"){

        double d = 2.116115162;
        rowvec posN1 = {-0.5*d, 0.0, 0.0};
        rowvec posN2= {0.5*d, 0.0, 0.0};
        rowvec charges = {7.0, 7.0};
        int nElectrons = 14;

        mat nucleiPositions = zeros<mat>(2,3);
        nucleiPositions.row(0) = posN1;
        nucleiPositions.row(1) = posN2;

        BasisFunctions* basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/N_6311Gs.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/N_6311Gs.dat", 1);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        RMP solver(system,2);
        solver.solve();

        if (my_rank == 0){
            cout << "Energy: " << setprecision(9) << solver.getEnergyHF() + solver.getEnergy2order() + solver.getEnergy3order() << endl;
        }

        delete system;
        delete basisFunctions;

    } else if (run == "FCl") {

        double d = 3.154519593;
        rowvec posF = {-0.5*d, 0.0, 0.0};
        rowvec posCl= {0.5*d, 0.0, 0.0};
        rowvec charges = {9.0, 17.0};
        int nElectrons = 26;

        mat nucleiPositions = zeros<mat>(2,3);
        nucleiPositions.row(0) = posF;
        nucleiPositions.row(1) = posCl;

        BasisFunctions* basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/F_6311Gs.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/Cl_6311Gs.dat", 1);

        System *system;
        system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

        RMP solver(system,2);
        solver.solve();

        if (my_rank == 0){
            cout << "Energy: " << setprecision(9) << solver.getEnergyHF() + solver.getEnergy2order() << endl;
        }

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
        BasisFunctions* basisFunctions = new BasisFunctions;
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/O_431G.dat", 0);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_431G.dat", 1);
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/H_431G.dat", 2);
        System *system = new System(basisFunctions, nucleiPositions, charges, nElectrons);
        RMP *solver = new RMP(system,1);
        HartreeFockFunc *func = new HartreeFockFunc(solver, system);
        Minimizer *minimizer = new Minimizer(func);
        minimizer->solve();

        if (my_rank == 0){
            cout << system->getNucleiPositions() << endl;
            cout << setprecision(7) << system->getNucleiPositions()(1,0) << endl;
            cout << "Energy min.: " << setprecision(14) << minimizer->getMinValue() << endl;
            cout << "Energy max.: " << setprecision(14) << minimizer->getMaxValue() << endl;
        }

    } else {
        if (my_rank == 0){
            cout << "No valid run selected." << endl;
        }
    }

    clock_t end = clock();
    if (my_rank == 0){
        cout << "Elapsed time: "<< (double(end - begin))/CLOCKS_PER_SEC << endl;
    }

#ifdef RUN_MPI
    MPI_Finalize();
#endif

    return 0;
}

