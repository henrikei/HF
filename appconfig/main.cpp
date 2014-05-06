#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <yaml-cpp/yaml.h>
#ifdef RUN_MPI
#include <mpi.h>
#endif
#include "hartreefock/rhf.h"
#include "hartreefock/uhf.h"
#include "perturbation/rmp.h"
#include "perturbation/ump.h"
#include "system/system.h"
#include "basisfunctions/basisfunctions.h"
#include "minimizer/minimizer.h"
#include "minimizer/hartreefockfunc.h"
#include "density/density.h"

using namespace std;
using namespace arma;
using namespace YAML;



struct Atom{
    string name;
    double charge;
    rowvec3 position;
};

void operator >> (const Node& node, rowvec3& position);
void operator >> (const Node& node, Atom& atom);

System* setup_system(const Node& doc);

void run_single(const Node& doc, char* argv[]);
//void run_multiple(const Node& doc, char* argv[]);
void run_minimize(const Node& doc, char* argv[]);

void write_density(const Node& node, field<mat> P, BasisFunctions *basisFunctions, char* argv[]);




int main(int argc, char* argv[])
{
#ifdef RUN_MPI
    MPI_Init(NULL, NULL);
#endif

    ifstream file;
    file.open(argv[1]);
    Parser parser(file);
    Node doc;
    parser.GetNextDocument(doc);

    string run_type;
    doc["run_type"] >> run_type;

    if (run_type == "single"){
        run_single(doc, argv);
    } else if (run_type == "minimize"){
        run_minimize(doc, argv);
    }

#ifdef RUN_MPI
    MPI_Finalize();
#endif

    return 0;
}




//---------------------------------------------------------------------------------------------------------------
void operator >> (const Node& node, rowvec3& position){
    node[0] >> position(0);
    node[1] >> position(1);
    node[2] >> position(2);
}


//-----------------------------------------------------------------------------------------------------------------
void operator >> (const Node& node, Atom& atom){
    node["name"] >> atom.name;
    node["charge"] >> atom.charge;
    node["position"] >> atom.position;
}


//------------------------------------------------------------------------------------------------------------------
System* setup_system(const Node& doc){
    vector<Atom> atoms;

    const Node& node = doc["atoms"];
    int nAtoms = node.size();
    for (int i = 0; i < nAtoms; i++){
        Atom atom;
        node[i] >> atom;
        atoms.push_back(atom);
    }

    mat nucleiPositions = zeros(nAtoms, 3);
    for (int i = 0; i < nAtoms; i++){
        nucleiPositions.row(i) = atoms[i].position;
    }

    BasisFunctions *basisFunctions = new BasisFunctions;
    string basis;
    doc["basis"] >> basis;
    for (int i = 0; i < nAtoms; i++){
        basisFunctions->addContracteds("../../HartreeFock/inFiles/basisSets/"+atoms[i].name+"_"+basis+".dat", i);
    }

    rowvec charges(nAtoms);
    for (int i = 0; i < nAtoms; i++){
        charges[i] = atoms[i].charge;
    }

    int n_electrons;
    doc["n_electrons"] >> n_electrons;

    System *system = new System(basisFunctions, nucleiPositions, charges, n_electrons);
    return system;
}


//---------------------------------------------------------------------------------------------------------------------
void run_single(const Node& doc, char *argv[]){

    System *system = setup_system(doc);
    ofstream file;
    string file_name = argv[2];
    file_name = file_name+"/energy.dat";
    file.open(file_name);

    int my_rank = 0;
#ifdef RUN_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

    string solver_type;
    int perturbation_order;
    doc["solver_type"] >> solver_type;
    if (solver_type == "RHF"){
        RHF solver(system);
        solver.solve();

        if (my_rank == 0){
            file << "Energy: " << setprecision(10) << solver.getEnergy() << endl;
            BasisFunctions *basisFunctions = system->getBasisFunctions();
            if (doc.FindValue("density")){
                const Node &node = doc["density"];
                field<mat> P = solver.getDensityMatrix();
                write_density(node, P, basisFunctions, argv);
            }
        }
    } else if (solver_type == "UHF"){
        UHF solver(system);
        solver.solve();

        if (my_rank == 0){
            file << "Energy: " << setprecision(10) << solver.getEnergy() << endl;
            BasisFunctions *basisFunctions = system->getBasisFunctions();
            if (doc.FindValue("density")){
                const Node &node = doc["density"];
                field<mat> P = solver.getDensityMatrix();
                write_density(node, P, basisFunctions, argv);
            }
        }
    } else if (solver_type == "RMP"){
        doc["perturbation_order"] >> perturbation_order;
        RMP solver(system, perturbation_order);
        solver.solve();

        if (my_rank == 0){
            file << "RHF energy: " << setprecision(10) << solver.getEnergyHF() << endl;
            file << "RMP2 energy: " << setprecision(10) << solver.getEnergyHF() + solver.getEnergy2order() << endl;
            file << "RMP3 energy: " << setprecision(10) << solver.getEnergyHF() + solver.getEnergy2order() + solver.getEnergy3order() << endl;
        }
    } else if (solver_type == "UMP"){
        doc["perturbation_order"] >> perturbation_order;
        UMP solver(system, perturbation_order);
        solver.solve();

        if (my_rank == 0){
            file << "UHF energy: " << setprecision(10) << solver.getEnergyHF() << endl;
            file << "UMP2 energy: " << setprecision(10) << solver.getEnergyHF() + solver.getEnergy2order() << endl;
            file << "UMP3 energy: " << setprecision(10) << solver.getEnergyHF() + solver.getEnergy2order() + solver.getEnergy3order() << endl;
        }
    } else {
        cout << "Error: Unknown solver type." << endl;
        exit(EXIT_FAILURE);
    }
    file.close();
}


//--------------------------------------------------------------------------------------------------------------
//void run_multiple(const Node& doc, char *argv[]){

//}


//---------------------------------------------------------------------------------------------------------------
void run_minimize(const Node& doc, char *argv[]){

    System *system = setup_system(doc);
    ofstream file;
    string file_name = argv[2];
    file_name = file_name+"/min_config.dat";
    file.open(file_name);

    string solver_type;
    int perturbation_order;
    doc["solver_type"] >> solver_type;

    int my_rank = 0;
#ifdef RUN_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif

    if (solver_type == "RHF"){

        RMP *solver = new RMP(system,1);
        HartreeFockFunc *func = new HartreeFockFunc(solver, system);
        Minimizer minimizer(func);
        minimizer.solve();

        if (my_rank == 0){
            file << "positions: " << endl;
            mat pos = func->getNucleiPositions();
            for (uint i = 0; i < pos.n_rows; i++){
                for (uint j = 0; j < pos.n_cols; j++){
                    file << setprecision(10) << pos(i,j) <<"  ";
                }
                file << endl;
            }
            file << endl << "Energy: " << setprecision(10) << minimizer.getMinValue() << endl;

            BasisFunctions *basisFunctions = system->getBasisFunctions();
            if (doc.FindValue("density")){
                const Node &node = doc["density"];
                field<mat> P = solver->getDensityMatrix();
                write_density(node, P, basisFunctions, argv);
            }
        }

    } else if (solver_type == "UHF"){

        UMP *solver = new UMP(system,1);
        HartreeFockFunc *func = new HartreeFockFunc(solver, system);
        Minimizer minimizer(func);
        minimizer.solve();

        if (my_rank == 0){
            file << "positions: " << endl;
            mat pos = func->getNucleiPositions();
            for (uint i = 0; i < pos.n_rows; i++){
                for (uint j = 0; j < pos.n_cols; j++){
                    file << setprecision(10) << pos(i,j) <<"  ";
                }
                file << endl;
            }
            file << endl  << "Energy: " << setprecision(10) << minimizer.getMinValue() << endl;

            BasisFunctions *basisFunctions = system->getBasisFunctions();
            if (doc.FindValue("density")){
                const Node &node = doc["density"];
                field<mat> P = solver->getDensityMatrix();
                write_density(node, P, basisFunctions, argv);
            }
        }

    } else if (solver_type == "RMP"){

        doc["perturbation_order"] >> perturbation_order;
        RMP *solver = new RMP(system, perturbation_order);
        HartreeFockFunc *func = new HartreeFockFunc(solver, system);
        Minimizer minimizer(func);
        minimizer.solve();

        if (my_rank == 0){
            file << "positions: " << endl;
            mat pos = func->getNucleiPositions();
            for (uint i = 0; i < pos.n_rows; i++){
                for (uint j = 0; j < pos.n_cols; j++){
                    file << setprecision(10) << pos(i,j) <<"  ";
                }
                file << endl;
            }
            file << endl << "Energy: " << setprecision(10) << minimizer.getMinValue() << endl;
        }

    } else if (solver_type == "UMP"){
        doc["perturbation_order"] >> perturbation_order;
        UMP *solver = new UMP(system, perturbation_order);
        HartreeFockFunc *func = new HartreeFockFunc(solver, system);
        Minimizer minimizer(func);
        minimizer.solve();

        if (my_rank == 0){
            file << "positions: " << endl;
            mat pos = func->getNucleiPositions();
            for (uint i = 0; i < pos.n_rows; i++){
                for (uint j = 0; j < pos.n_cols; j++){
                    file << setprecision(10) << pos(i,j) <<"  ";
                }
                file << endl;
            }
            file << endl << "Energy: " << setprecision(10) << minimizer.getMinValue() << endl;
        }

    } else {
        cout << "Error: Unknown solver type." << endl;
        exit(EXIT_FAILURE);
    }
    file.close();
}


//---------------------------------------------------------------------------------------------------------------------
void write_density(const Node& node, field<mat> P, BasisFunctions *basisFunctions, char* argv[]){
    rowvec3 R1, R2;
    double dx, dy, dz;
    node["corner1"] >> R1;
    node["corner2"] >> R2;
    node["dx"] >> dx; node["dy"] >> dy; node["dz"] >> dz;

    Density density(basisFunctions, R1, R2, dx, dy, dz);
    string filename = argv[2];
    filename = filename+"/density.dat";
    density.printDensity(P, filename);
}
