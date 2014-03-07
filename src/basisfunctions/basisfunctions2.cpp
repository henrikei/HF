#include "basisfunctions2.h"

BasisFunctions2::BasisFunctions2()
{
}

//---------------------------------------------------------------------------------------------------------------------------
// Read exponents, coefficients and powers from input file in format TurboMole (taken from https://bse.pnl.gov/bse/portal)
// and fills m_contracteds (vector<Contracted*> m_contracteds) with contracteds. Each contracted contains a vector of primitive.
// Each primitive contains exponent, coefficient, powers (i,j,k) and position.
void BasisFunctions2::addContracteds(string inFileName, rowvec3 position)
{
    // pows_s, pows_p and pows_d defines the powers which are used to form primitives. For example, i+j+k = 1 for p-orbitals.
    // The prefactors are needed in cases where each primitive has a sum of two or more power terms, which is the case when
    // 5 d-orbitals are used instead of 6. The 5 p-orbitals are have the following powers:
    // xy, xz, yz, x^2 - y^2, 3*r^2 - z^2.

    field <imat> pows_s(1);
    field <imat> pows_p(3);
    field <imat> pows_d;
    field <rowvec> prefactor_s(1);
    field <rowvec> prefactor_p(3);
    field <rowvec> prefactor_d;

    bool d6 = true;
    if (d6 == true){
        pows_s(0) = {0,0,0};
        pows_p(0) = {1,0,0};
        pows_p(1) = {0,1,0};
        pows_p(2) = {0,0,1};
        pows_d.set_size(6);
        pows_d(0) = {2,0,0};
        pows_d(1) = {0,2,0};
        pows_d(2) = {0,0,2};
        pows_d(3) = {1,1,0};
        pows_d(4) = {1,0,1};
        pows_d(5) = {0,1,1};
        prefactor_s(0) = {1.0};
        prefactor_p(0) = {1.0};
        prefactor_p(1) = {1.0};
        prefactor_p(2) = {1.0};
        prefactor_d.set_size(6);
        prefactor_d(0) = {1.0};
        prefactor_d(1) = {1.0};
        prefactor_d(2) = {1.0};
        prefactor_d(3) = {1.0};
        prefactor_d(4) = {1.0};
        prefactor_d(5) = {1.0};
    } else {
        pows_s(0) = {0,0,0};
        pows_p(0) = {1,0,0};
        pows_p(1) = {0,1,0};
        pows_p(2) = {0,0,1};
        pows_d.set_size(5);
        pows_d(0) = {1,1,0};
        pows_d(1) = {1,0,1};
        pows_d(2) = {0,1,1};
        imat pows_dTemp = zeros<imat>(3,3);
        pows_dTemp.diag() += 2;
        pows_d(3) = pows_dTemp.rows(0,1);
        pows_d(4) = pows_dTemp;
        prefactor_s(0) = {1.0};
        prefactor_p(0) = {1.0};
        prefactor_p(1) = {1.0};
        prefactor_p(2) = {1.0};
        prefactor_d.set_size(5);
        prefactor_d(0) = {1.0};
        prefactor_d(1) = {1.0};
        prefactor_d(2) = {1.0};
        prefactor_d(3) = {1.0, -1.0};
        prefactor_d(4) = {-1.0, -1.0, 2.0};
    }


    string line, stringToSearch;
    fstream file;
    file.open(inFileName);
    if (file.is_open()){
        while (getline(file, line)){
            stringToSearch += line +"\n";
        }
        file.close();
    } else {
        cout << "Error: Could not open file "<< inFileName << endl;
        exit(EXIT_FAILURE);
    }

    // search1: Search for groups of contracted basis sets
    regex search1("([0-9]\\s+[spd])((\\s+-?[0-9]+\\.[0-9]+)+)");
    sregex_iterator pos1(stringToSearch.begin(), stringToSearch.end(), search1);
    sregex_iterator end1;
    for(; pos1!=end1; pos1++){
        string subString = pos1->str(2).c_str();
        // search2: Search for groups of (exponent and coefficient) for each primitive
        regex search2("(-?[0-9]+\\.[0-9]+)\\s+(-?[0-9]+\\.[0-9]+)\\s*");
        sregex_iterator pos2(subString.begin(), subString.end(), search2);
        sregex_iterator end2;
        vector<double> exp;
        vector<double> coeff;
        // collect exponent and coefficient for each primitive
        for(; pos2!=end2; pos2++){
            exp.push_back(atof(pos2->str(1).c_str()));
            coeff.push_back(atof(pos2->str(2).c_str()));
        }

        // Check type of orbital (s, p, d)
        string orbitalType = pos1->str(1).c_str();
        regex search_s("[0-9]\\s+s");
        regex search_p("[0-9]\\s+p");
        regex search_d("[0-9]\\s+d");
        if(regex_match(orbitalType, search_s)){
            addSomeContracteds(exp, coeff, pows_s, prefactor_s, position);
        }
        if(regex_match(orbitalType, search_p)){
            addSomeContracteds(exp, coeff, pows_p, prefactor_p, position);
        }
        if(regex_match(orbitalType, search_d)){
            addSomeContracteds(exp, coeff, pows_d, prefactor_d, position);
        }
    }
    //----------------------------------------------------------------------------------------------------------------------------
}

double BasisFunctions2::getNumOfContracteds()
{
    return m_contracteds.size();
}

Contracted *BasisFunctions2::getContracted(int p)
{
    return m_contracteds.at(p);
}

int BasisFunctions2::getAngMomMax()
{
    int maxAngMom = 0;
    for (uint i = 0; i < m_contracteds.size(); i++){
        if (maxAngMom < m_contracteds.at(i)->getAngMom()){
            maxAngMom = m_contracteds.at(i)->getAngMom();
        }
    }
    return maxAngMom;
}

void BasisFunctions2::addSomeContracteds(vector<double> exp, vector<double> coeff, field<imat> pows, field<rowvec> prefactor, rowvec3 position)
{
    // Loop through number of contracteds (1s -> 1, 2p -> 3, 2d -> 6 (or 5)
    for (uint i = 0; i < pows.n_rows; i++){
        vector<Primitive*> primitives;
        // Loop though number of power terms (only different from 1 when 5 d-functions are used)
        for (uint j = 0; j < pows(i).n_rows; j++){
            int counter = 0;
            // Loop through the elements of exp and coeff
            while (counter < (int)exp.size()){
                double coeffn = coeff.at(counter);
                normalizeCoeff(exp.at(counter), coeffn, pows(i).row(j));
                coeffn *= prefactor(i)(j);
                Primitive *primitive = new Primitive(exp.at(counter), coeffn, pows(i).row(j), position);
                primitives.push_back(primitive);
                counter += 1;
            }
        }
        Contracted *contracted = new Contracted(primitives);
        m_contracteds.push_back(contracted);
    }
}

void BasisFunctions2::normalizeCoeff(double exp, double &coeff, irowvec pows)
{
    int i = pows.at(0);
    int j = pows.at(1);
    int k = pows.at(2);
    coeff *= pow((2*exp/M_PI),0.75)*sqrt(pow(8*exp,i+j+k)*factorial(i)*factorial(j)*factorial(k)/(factorial(2*i)*factorial(2*j)*factorial(2*k)));
}

int BasisFunctions2::factorial(int n)
{
    double value = 1;
    double i = 1;

    while(i < n){
        i += 1;
        value *= i;
    }
    return value;
}
