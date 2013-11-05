#include "o_321g.h"

O_321G::O_321G()
{
    rowvec exp1 = {322.0370000, 48.430800, 10.4206000};
    rowvec exp2 = {7.4029400, 1.5762000};
    rowvec exp3 = {0.3736840};
    rowvec exp4 = {7.4029400, 1.5762000};
    rowvec exp5 = {0.3736840};

    rowvec coeffs1 = {0.0592394, 0.3515000, 0.7076580};
    rowvec coeffs2 = {-0.4044530, 1.2215600};
    rowvec coeffs3 = {1.0};
    rowvec coeffs4 = {0.2445860, 0.8539550};
    rowvec coeffs5 = {1.0};

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
}
