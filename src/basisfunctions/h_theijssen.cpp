#include "h_theijssen.h"

H_Theijssen::H_Theijssen()
{
    rowvec exp1 = {13.00773};
    rowvec exp2 = {1.962079};
    rowvec exp3 = {0.444529};
    rowvec exp4 = {0.1219492};

    rowvec coeffs1 = {1.0};
    rowvec coeffs2 = {1.0};
    rowvec coeffs3 = {1.0};
    rowvec coeffs4 = {1.0};

    irowvec powers1 = {0, 0, 0};
    irowvec powers2 = {0, 0, 0};
    irowvec powers3 = {0, 0, 0};
    irowvec powers4 = {0, 0, 0};

    exponents.push_back(exp1);
    exponents.push_back(exp2);
    exponents.push_back(exp3);
    exponents.push_back(exp4);

    coeffs.push_back(coeffs1);
    coeffs.push_back(coeffs2);
    coeffs.push_back(coeffs3);
    coeffs.push_back(coeffs4);

    powers.push_back(powers1);
    powers.push_back(powers2);
    powers.push_back(powers3);
    powers.push_back(powers4);

    angMom = 0;
}

