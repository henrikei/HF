#include "basisfunctions.h"


BasisFunctions::BasisFunctions()
{
}

BasisFunctions::BasisFunctions(string inFileName)
{
    //---------------------------------------------------------------------------------------------------------------------------
    // Read exponents, coefficients and powers from input file in format TurboMole (taken from https://bse.pnl.gov/bse/portal)
    irowvec pow1 = {0, 0, 0};
    irowvec pow2 = {1, 0, 0};
    irowvec pow3 = {0, 1, 0};
    irowvec pow4 = {0, 0, 1};
    irowvec pow5 = {2, 0, 0};
    irowvec pow6 = {0, 2, 0};
    irowvec pow7 = {0, 0, 2};
    irowvec pow8 = {1, 1, 0};
    irowvec pow9 = {1, 0, 1};
    irowvec pow10 = {0, 1, 1};

    string line, stringToSearch;
    fstream file;
    file.open(inFileName);
    if (file.is_open()){
        while (getline(file, line)){
            stringToSearch += line +"\n";
        }
        file.close();
    } else {
        cout << "Error: Could not open file." << endl;
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
        vector<double> temp1;
        vector<double> temp2;
        for(; pos2!=end2; pos2++){
            temp1.push_back(atof(pos2->str(1).c_str()));
            temp2.push_back(atof(pos2->str(2).c_str()));
        }
        rowvec exp = zeros<rowvec>(temp1.size());
        rowvec c = zeros<rowvec>(temp1.size());
        // Fill exp and c with exponents and coefficients of current basis
        for(uint i = 0; i < temp1.size(); i++){
            exp(i) = temp1.at(i);
            c(i) = temp2.at(i);
        }

        // Check type of orbital (s, p, d)
        string orbitalType = pos1->str(1).c_str();
        regex searchs("[0-9]\\s+s");
        regex searchp("[0-9]\\s+p");
        regex searchd("[0-9]\\s+d");
        if(regex_match(orbitalType, searchs)){
            exponents.push_back(exp);
            coeffs.push_back(c);
            powers.push_back(pow1);
        }
        if(regex_match(orbitalType, searchp)){
            exponents.push_back(exp);
            exponents.push_back(exp);
            exponents.push_back(exp);
            coeffs.push_back(c);
            coeffs.push_back(c);
            coeffs.push_back(c);
            powers.push_back(pow2);
            powers.push_back(pow3);
            powers.push_back(pow4);
        }
        if(regex_match(orbitalType, searchd)){
            exponents.push_back(exp);
            exponents.push_back(exp);
            exponents.push_back(exp);
            exponents.push_back(exp);
            exponents.push_back(exp);
            exponents.push_back(exp);
            coeffs.push_back(c);
            coeffs.push_back(c);
            coeffs.push_back(c);
            coeffs.push_back(c);
            coeffs.push_back(c);
            coeffs.push_back(c);
            powers.push_back(pow5);
            powers.push_back(pow6);
            powers.push_back(pow7);
            powers.push_back(pow8);
            powers.push_back(pow9);
            powers.push_back(pow10);
        }
    }
    //----------------------------------------------------------------------------------------------------------------------------

    normalizeCoeffs();

    // Find max angular momentum of basis set.
    angMom = 0;
    for (uint i = 0; i < powers.size(); i++){
        if(angMom < sum(powers.at(i))){
            angMom = sum(powers.at(i));
        }
    }
}


void BasisFunctions::setPosition(rowvec newPosition)
{
    position = newPosition;
}

rowvec3 BasisFunctions::getPosition()
{
    return position;
}


rowvec BasisFunctions::getExponents(int p)
{
    return exponents.at(p);
}


rowvec BasisFunctions::getCoeffs(int p)
{
    return coeffs.at(p);
}

irowvec BasisFunctions::getPowers(int p)
{
    return powers.at(p);
}

int BasisFunctions::getAngMom(){
    return angMom;
}

int BasisFunctions::getNumOfBasisFunc()
{
    return exponents.size();
}

void BasisFunctions::normalizeCoeffs(){
    for (int i = 0; i < (int)coeffs.size(); i++){
        int l = sum(powers.at(i));
        for (int j = 0; j < (int)coeffs.at(i).n_elem; j++){
            if (l == 0){
                coeffs.at(i)(j) *= pow(2*exponents.at(i)(j)/M_PI, 0.75);
            } else if (l == 1){
                coeffs.at(i)(j) *= pow(2*exponents.at(i)(j)/M_PI, 0.75)*2*sqrt(exponents.at(i)(j));
            } else if (l == 2){
                int lx = powers.at(i)(0);
                int ly = powers.at(i)(1);
                int lz = powers.at(i)(2);
                coeffs.at(i)(j) *= pow(2*exponents.at(i)(j)/M_PI, 0.75)*sqrt(pow(8*exponents.at(i)(j),lx+ly+lz)
                                   *factorial(2*lx)*factorial(2*ly)*factorial(2*lz)/(factorial(2*lx)*factorial(2*ly)*factorial(2*lz)));
            }
        }
    }
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
