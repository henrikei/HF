#include "h_6311gss.h"

H_6311Gss::H_6311Gss()
{
    rowvec exp1 = {33.8650000, 5.0947900, 1.1587900};
    rowvec exp2 = {0.3258400};
    rowvec exp3 = {0.1027410};
    rowvec exp4 = {0.7500000};
    rowvec exp5 = {0.7500000};
    rowvec exp6 = {0.7500000};

    rowvec coeffs1 = {0.0254938, 0.1903730, 0.8521610};
    rowvec coeffs2 = {1.0};
    rowvec coeffs3 = {1.0};
    rowvec coeffs4 = {1.0};
    rowvec coeffs5 = {1.0};
    rowvec coeffs6 = {1.0};

    irowvec powers1 = {0, 0, 0};
    irowvec powers2 = {0, 0, 0};
    irowvec powers3 = {0, 0, 0};
    irowvec powers4 = {1, 0, 0};
    irowvec powers5 = {0, 1, 0};
    irowvec powers6 = {0, 0, 1};


    exponents.push_back(exp1);
    exponents.push_back(exp2);
    exponents.push_back(exp3);
    exponents.push_back(exp4);
    exponents.push_back(exp5);
    exponents.push_back(exp6);

    coeffs.push_back(coeffs1);
    coeffs.push_back(coeffs2);
    coeffs.push_back(coeffs3);
    coeffs.push_back(coeffs4);
    coeffs.push_back(coeffs5);
    coeffs.push_back(coeffs6);

    powers.push_back(powers1);
    powers.push_back(powers2);
    powers.push_back(powers3);
    powers.push_back(powers4);
    powers.push_back(powers5);
    powers.push_back(powers6);

    angMom = 1;

    normalizeCoeffs();
}
