#include <iostream>
#include <iomanip>
#include <fstream>
#include <armadillo>
#include <yaml-cpp/yaml.h>
#include "solverwrapper/solverwrapper.h"
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

void run_single(const Node& doc);
void run_multiple(const Node& doc);
void run_minimize(const Node& doc);




int main()
{
    ifstream file;
    file.open("../../HartreeFock/appconfig/configs/test.cfg");
    Parser parser(file);
    Node doc;
    parser.GetNextDocument(doc);

    string run_type;
    doc["run_type"] >> run_type;

    if (run_type == "single"){
        run_single(doc);
    } else if (run_type == "multiple"){
        run_multiple(doc);
    } else if (run_type == "minimize"){
        run_minimize(doc);
    }

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
void run_single(const Node& doc){

    System *system = setup_system(doc);
    ofstream file;
    file.open("../../Results/test.out");

    string solver_type;
    int perturbation_order;
    doc["solver_type"] >> solver_type;
    if (solver_type == "RHF"){
        SolverWrapper<RHF> solver(system);
        solver.solve();
        file << "Energy: " << solver.getEnergy() << endl;
    } else if (solver_type == "UHF"){
        SolverWrapper<UHF> solver(system);
        solver.solve();
        file << "Energy: " << solver.getEnergy() << endl;
    } else if (solver_type == "RMP"){
        doc["perturbation_order"] >> perturbation_order;
        SolverWrapper<RMP> solver(system, perturbation_order);
        solver.solve();
        file << "Energy: " << solver.getEnergy() << endl;
    } else if (solver_type == "UMP"){
        doc["perturbation_order"] >> perturbation_order;
        SolverWrapper<UMP> solver(system, perturbation_order);
        solver.solve();
        file << "Energy: " << solver.getEnergy() << endl;
    } else {
        cout << "Error: Unknown solver type." << endl;
        exit(EXIT_FAILURE);
    }
    file.close();
}


//--------------------------------------------------------------------------------------------------------------
void run_multiple(const Node& doc){

}


//---------------------------------------------------------------------------------------------------------------
void run_minimize(const Node& doc){

    System *system = setup_system(doc);
    ofstream file;
    file.open("../../Results/test.out");

    string solver_type;
    int perturbation_order;
    doc["solver_type"] >> solver_type;
    if (solver_type == "RHF"){
        RMP *solver = new RMP(system,1);
        HartreeFockFunc *func = new HartreeFockFunc(solver, system);
        Minimizer minimizer(func);
        minimizer.solve();
        mat pos = system->getNucleiPositions();
        for (uint i = 0; i < pos.n_rows; i++){
            for (uint j = 0; j < pos.n_cols; j++){
                file << setprecision(10) << pos(i,j) <<"  ";
            }
            file << endl;
        }
        file << endl << minimizer.getMinValue() << endl;
    } else if (solver_type == "UHF"){
        UMP *solver = new UMP(system,1);
        HartreeFockFunc *func = new HartreeFockFunc(solver, system);
        Minimizer minimizer(func);
        minimizer.solve();
        mat pos = system->getNucleiPositions();
        for (int i = 0; i < pos.n_rows; i++){
            for (int j = 0; j < pos.n_cols; j++){
                file << setprecision(10) << pos(i,j) <<"  ";
            }
            file << endl;
        }
        file << endl << minimizer.getMinValue() << endl;
    } else if (solver_type == "RMP"){
        doc["perturbation_order"] >> perturbation_order;
        RMP *solver = new RMP(system, perturbation_order);
        HartreeFockFunc *func = new HartreeFockFunc(solver, system);
        Minimizer minimizer(func);
        minimizer.solve();
        mat pos = system->getNucleiPositions();
        for (int i = 0; i < pos.n_rows; i++){
            for (int j = 0; j < pos.n_cols; j++){
                file << setprecision(10) << pos(i,j) <<"  ";
            }
            file << endl;
        }
        file << endl << minimizer.getMinValue() << endl;
    } else if (solver_type == "UMP"){
        doc["perturbation_order"] >> perturbation_order;
        UMP *solver = new UMP(system, perturbation_order);
        HartreeFockFunc *func = new HartreeFockFunc(solver, system);
        Minimizer minimizer(func);
        minimizer.solve();
        mat pos = system->getNucleiPositions();
        for (int i = 0; i < pos.n_rows; i++){
            for (int j = 0; j < pos.n_cols; j++){
                file << setprecision(10) << pos(i,j) <<"  ";
            }
            file << endl;
        }
        file << endl << minimizer.getMinValue() << endl;
    } else {
        cout << "Error: Unknown solver type." << endl;
        exit(EXIT_FAILURE);
    }
    file.close();
}
