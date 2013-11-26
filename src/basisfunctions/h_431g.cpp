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

    normalizeCoeffs();
}
