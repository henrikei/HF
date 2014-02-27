#include "basisfunctions2.h"

BasisFunctions2::BasisFunctions2()
{
}

void BasisFunctions2::addContracteds(string inFileName, rowvec3 position)
{
    //---------------------------------------------------------------------------------------------------------------------------
    // Read exponents, coefficients and powers from input file in format TurboMole (taken from https://bse.pnl.gov/bse/portal)
    imat pows = zeros<imat>(10,3);
    pows(0,0) = 0; pows(0,1) = 0; pows(0,2) = 0;
    pows(1,0) = 1; pows(1,1) = 0; pows(1,2) = 0;
    pows(2,0) = 0; pows(2,1) = 1; pows(2,2) = 0;
    pows(3,0) = 0; pows(3,1) = 0; pows(3,2) = 1;
    pows(4,0) = 2; pows(4,1) = 0; pows(4,2) = 0;
    pows(5,0) = 0; pows(5,1) = 2; pows(5,2) = 0;
    pows(6,0) = 0; pows(6,1) = 0; pows(6,2) = 2;
    pows(7,0) = 1; pows(7,1) = 1; pows(7,2) = 0;
    pows(8,0) = 1; pows(8,1) = 0; pows(8,2) = 1;
    pows(9,0) = 0; pows(9,1) = 1; pows(9,2) = 1;

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
        regex searchs("[0-9]\\s+s");
        regex searchp("[0-9]\\s+p");
        regex searchd("[0-9]\\s+d");
        if(regex_match(orbitalType, searchs)){
            int counter = 0;
            vector<Primitive*> primitives;
            while (counter < (int)exp.size()){
                double coeffn = coeff.at(counter);
                normalizeCoeff(exp.at(counter), coeffn, pows.row(0));
                Primitive *primitive = new Primitive(exp.at(counter), coeffn, pows.row(0), position);
                primitives.push_back(primitive);
                counter += 1;
            }
            Contracted *contracted = new Contracted(primitives);
            m_contracteds.push_back(contracted);
        }
        if(regex_match(orbitalType, searchp)){
            for (int i = 1; i < 4; i++){
                int counter = 0;
                vector<Primitive*> primitives;
                while (counter < (int)exp.size()){
                    double coeffn = coeff.at(counter);
                    normalizeCoeff(exp.at(counter), coeffn, pows.row(i));
                    Primitive *primitive = new Primitive(exp.at(counter), coeffn, pows.row(i), position);
                    primitives.push_back(primitive);
                    counter += 1;
                }
                Contracted *contracted = new Contracted(primitives);
                m_contracteds.push_back(contracted);
            }
        }
        if(regex_match(orbitalType, searchd)){
            for (int i = 4; i < 10; i++){
                int counter = 0;
                vector<Primitive*> primitives;
                while (counter < (int)exp.size()){
                    double coeffn = coeff.at(counter);
                    normalizeCoeff(exp.at(counter), coeffn, pows.row(i));
                    Primitive *primitive = new Primitive(exp.at(counter), coeffn, pows.row(i), position);
                    primitives.push_back(primitive);
                    counter += 1;
                }
                Contracted *contracted = new Contracted(primitives);
                m_contracteds.push_back(contracted);
            }
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
    for (int i = 0; i < (int)m_contracteds.size(); i++){
        if (maxAngMom < m_contracteds.at(i)->getAngMom()){
            maxAngMom = m_contracteds.at(i)->getAngMom();
        }
    }
    return maxAngMom;
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
