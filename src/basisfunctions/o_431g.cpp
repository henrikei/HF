#include "o_431g.h"

O_431G::O_431G()
{
    rowvec exp1 = {883.2728600, 133.1292800, 29.9064080, 7.9786772};
    rowvec exp2 = {16.1944470, 3.7800860, 1.0709836};
    rowvec exp3 = {0.2838798};
    rowvec exp4 = {16.1944470, 3.7800860, 1.0709836};
    rowvec exp5 = {0.2838798};

    rowvec coeffs1 = {0.0175506, 0.1228292, 0.4348836, 0.5600108};
    rowvec coeffs2 = {-0.1134010, -0.1772865, 1.1504079};
    rowvec coeffs3 = {1.0000000};
    rowvec coeffs4 = {0.0685453, 0.3312254 ,0.7346079};
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

    normalizeCoeffs();
}
