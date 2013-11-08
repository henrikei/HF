#include "n_431g.h"

N_431G::N_431G()
{
    rowvec exp1 = {671.2795000, 101.2017000, 22.6999700, 6.0406090};
    rowvec exp2 = {12.3935997, 2.9223828, 0.83252808};
    rowvec exp3 = {0.2259640};
    rowvec exp4 = {12.3935997, 2.9223828, 0.83252808};
    rowvec exp5 = {0.2259640};

    rowvec coeffs1 = {0.0175982511, 0.1228462410, 0.4337821410, 0.5614182170};
    rowvec coeffs2 = {-0.1174892990, -0.2139940160, 1.1745021100};
    rowvec coeffs3 = {1.0000000};
    rowvec coeffs4 = {0.0640203443, 0.3112025550, 0.7527482390};
    rowvec coeffs5 = {1.0000000};

    irowvec powers1 = {0, 0, 0};
    irowvec powers2 = {1, 0, 0};
    irowvec powers3 = {0, 1, 0};
    irowvec powers4 = {0, 0, 1};

    exponents.push_back(exp1);
    exponents.push_back(exp2);
    exponents.push_back(exp3);
    exponents.push_back(exp4);
    exponents.push_back(exp5);
    exponents.push_back(exp4);
    exponents.push_back(exp5);
    exponents.push_back(exp4);
    exponents.push_back(exp5);

    coeffs.push_back(coeffs1);
    coeffs.push_back(coeffs2);
    coeffs.push_back(coeffs3);
    coeffs.push_back(coeffs4);
    coeffs.push_back(coeffs5);
    coeffs.push_back(coeffs4);
    coeffs.push_back(coeffs5);
    coeffs.push_back(coeffs4);
    coeffs.push_back(coeffs5);

    powers.push_back(powers1);
    powers.push_back(powers1);
    powers.push_back(powers1);
    powers.push_back(powers2);
    powers.push_back(powers2);
    powers.push_back(powers3);
    powers.push_back(powers3);
    powers.push_back(powers4);
    powers.push_back(powers4);

    angMom = 1;

    for (int i = 0; i < (int)coeffs.size(); i++){
        int l = sum(powers.at(i));
        for (int j = 0; j < (int)coeffs.at(i).n_elem; j++){
            if (l == 0){
                coeffs.at(i)(j) *= pow(2*exponents.at(i)(j)/M_PI, 0.75);
            } else if (l == 1){
                coeffs.at(i)(j) *= pow(2*exponents.at(i)(j)/M_PI, 0.75)*2*sqrt(exponents.at(i)(j));
            }
        }
    }
}
