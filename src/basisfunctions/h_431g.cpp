#include "h_431g.h"

H_431G::H_431G()
{
    rowvec exp1 = {18.7311370, 2.8253944, 0.6401217};
    rowvec exp2 = {0.1612778};

    rowvec coeffs1 = {0.0334946, 0.2347269, 0.8137573};
    rowvec coeffs2 = {1.0};

    irowvec powers1 = {0, 0, 0};
    irowvec powers2 = {0, 0, 0};


    exponents.push_back(exp1);
    exponents.push_back(exp2);

    coeffs.push_back(coeffs1);
    coeffs.push_back(coeffs2);

    powers.push_back(powers1);
    powers.push_back(powers2);

    angMom = 0;

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
