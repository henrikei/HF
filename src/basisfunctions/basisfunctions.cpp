#include "basisfunctions.h"

BasisFunctions::BasisFunctions()
{
}

BasisFunctions::~BasisFunctions()
{
    for (Contracted* contracted: m_contracteds){
        for (int i = 0; i < contracted->getNumOfPrimitives(); i++){
            Primitive* primitive = contracted->getPrimitive(i);
            delete primitive;
        }
        delete contracted;
    }
}

//---------------------------------------------------------------------------------------------------------------------------
// Read exponents, coefficients and powers from input file in format TurboMole (taken from https://bse.pnl.gov/bse/portal)
// and fills m_contracteds (vector<Contracted*> m_contracteds) with contracteds. Each contracted contains a vector of primitive.
// Each primitive contains exponent, coefficient, powers (i,j,k) and position number.
void BasisFunctions::addContracteds(string inFileName, int posNum)
{
    // pows_s, pows_p and pows_d defines the powers which are used to form primitives. For example, i+j+k = 1 for p-orbitals.

    field<irowvec3> pows_s(1,3);
    field<irowvec3> pows_p(3,3);
    field<irowvec3> pows_d(6,3);
    field<irowvec3> pows_f(10,3);

    pows_s(0) = {0,0,0};
    pows_p(0) = {1,0,0};
    pows_p(1) = {0,1,0};
    pows_p(2) = {0,0,1};
    pows_d(0) = {2,0,0};
    pows_d(1) = {0,2,0};
    pows_d(2) = {0,0,2};
    pows_d(3) = {1,1,0};
    pows_d(4) = {1,0,1};
    pows_d(5) = {0,1,1};
    pows_f(0) = {3,0,0};
    pows_f(1) = {0,3,0};
    pows_f(2) = {0,0,3};
    pows_f(3) = {2,1,0};
    pows_f(4) = {2,0,1};
    pows_f(5) = {1,2,0};
    pows_f(6) = {0,2,1};
    pows_f(7) = {1,0,2};
    pows_f(8) = {0,1,2};
    pows_f(9) = {1,1,1};

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
    regex search1("([0-9]\\s+[spdf])((\\s+-?[0-9]+\\.[0-9]+)+)");
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

        // Check type of orbital (s, p, d, f)
        string orbitalType = pos1->str(1).c_str();
        regex search_s("[0-9]\\s+s");
        regex search_p("[0-9]\\s+p");
        regex search_d("[0-9]\\s+d");
        regex search_f("[0-9]\\s+f");
        if(regex_match(orbitalType, search_s)){
            addSomeContracteds(exp, coeff, pows_s, posNum);
        } else if (regex_match(orbitalType, search_p)){
            addSomeContracteds(exp, coeff, pows_p, posNum);
        } else if (regex_match(orbitalType, search_d)){
            addSomeContracteds(exp, coeff, pows_d, posNum);
        } else if (regex_match(orbitalType, search_f)){
            addSomeContracteds(exp, coeff, pows_f, posNum);
        }
    }
    //----------------------------------------------------------------------------------------------------------------------------
}

void BasisFunctions::setPosPointer(mat *nucleiPositions)
{
    for (Contracted* contracted: m_contracteds){
        for (int i = 0; i < contracted->getNumOfPrimitives(); i++){
            Primitive* primitive = contracted->getPrimitive(i);
            primitive->setPosPointer(nucleiPositions);
        }
    }
}

double BasisFunctions::getNumOfContracteds()
{
    return m_contracteds.size();
}

Contracted *BasisFunctions::getContracted(int p)
{
    return m_contracteds.at(p);
}

int BasisFunctions::getAngMomMax()
{
    int maxAngMom = 0;
    for (uint i = 0; i < m_contracteds.size(); i++){
        if (maxAngMom < m_contracteds.at(i)->getAngMom()){
            maxAngMom = m_contracteds.at(i)->getAngMom();
        }
    }
    return maxAngMom;
}

void BasisFunctions::addSomeContracteds(vector<double> exp, vector<double> coeff, field<irowvec3> pows, int posNum)
{
    // Loop through number of contracteds (s -> 1, p -> 3, d -> 6, f -> 10)
    for (uint i = 0; i < pows.n_rows; i++){
        vector<Primitive*> primitives;
        int counter = 0;
        // Loop through the elements of exp and coeff
        while (counter < (int)exp.size()){
            double coeffn = coeff.at(counter);
            normalizeCoeff(exp.at(counter), coeffn, pows(i));
            Primitive *primitive = new Primitive(exp.at(counter), coeffn, pows(i), posNum);
            primitives.push_back(primitive);
            counter += 1;
        }
        Contracted *contracted = new Contracted(primitives);
        m_contracteds.push_back(contracted);
    }
}

void BasisFunctions::normalizeCoeff(double exp, double &coeff, irowvec3 pows)
{
    int i = pows.at(0);
    int j = pows.at(1);
    int k = pows.at(2);
    coeff *= pow((2*exp/M_PI),0.75)*sqrt(pow(8*exp,i+j+k)*factorial(i)*factorial(j)*factorial(k)/(factorial(2*i)*factorial(2*j)*factorial(2*k)));
}

int BasisFunctions::factorial(int n)
{
    double value = 1;
    double i = 1;

    while(i < n){
        i += 1;
        value *= i;
    }
    return value;
}
