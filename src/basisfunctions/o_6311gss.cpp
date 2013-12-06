#include "o_6311gss.h"

O_6311Gss::O_6311Gss()
{
    rowvec exp1 = {8588.5000000, 1297.2300000, 299.2960000, 87.3771000, 25.6789000, 3.7400400};
    rowvec exp2 = {42.1175000, 9.6283700, 2.8533200};
    rowvec exp3 = {0.9056610};
    rowvec exp4 = {0.2556110};
    rowvec exp5 = {42.1175000, 9.6283700, 2.8533200};
    rowvec exp6 = {0.9056610};
    rowvec exp7 = {0.2556110};
    rowvec exp8 = {1.2920000};

    rowvec coeffs1 = {0.00189515, 0.0143859, 0.0707320, 0.2400010, 0.5947970, 0.2808020};
    rowvec coeffs2 = {0.1138890, 0.9208110, -0.00327447};
    rowvec coeffs3 = {1.0};
    rowvec coeffs4 = {1.0};
    rowvec coeffs5 = {0.0365114, 0.2371530, 0.8197020};
    rowvec coeffs6 = {1.0};
    rowvec coeffs7 = {1.0};
    rowvec coeffs8 = {1.0};

    irowvec powers1 = {0, 0, 0};
    irowvec powers2 = {1, 0, 0};
    irowvec powers3 = {0, 1, 0};
    irowvec powers4 = {0, 0, 1};
    irowvec powers5 = {1, 1, 0};
    irowvec powers6 = {1, 0, 1};
    irowvec powers7 = {0, 1, 1};
    irowvec powers8 = {2, 0, 0};
    irowvec powers9 = {0, 2, 0};
    irowvec powers10 = {0, 0, 2};

    exponents.push_back(exp1);
    exponents.push_back(exp2);
    exponents.push_back(exp3);
    exponents.push_back(exp4);
    exponents.push_back(exp5);
    exponents.push_back(exp6);
    exponents.push_back(exp7);
    exponents.push_back(exp5);
    exponents.push_back(exp6);
    exponents.push_back(exp7);
    exponents.push_back(exp5);
    exponents.push_back(exp6);
    exponents.push_back(exp7);
    exponents.push_back(exp8);
    exponents.push_back(exp8);
    exponents.push_back(exp8);
    exponents.push_back(exp8);
    exponents.push_back(exp8);
    exponents.push_back(exp8);

    coeffs.push_back(coeffs1);
    coeffs.push_back(coeffs2);
    coeffs.push_back(coeffs3);
    coeffs.push_back(coeffs4);
    coeffs.push_back(coeffs5);
    coeffs.push_back(coeffs6);
    coeffs.push_back(coeffs7);
    coeffs.push_back(coeffs5);
    coeffs.push_back(coeffs6);
    coeffs.push_back(coeffs7);
    coeffs.push_back(coeffs5);
    coeffs.push_back(coeffs6);
    coeffs.push_back(coeffs7);
    coeffs.push_back(coeffs8);
    coeffs.push_back(coeffs8);
    coeffs.push_back(coeffs8);
    coeffs.push_back(coeffs8);
    coeffs.push_back(coeffs8);
    coeffs.push_back(coeffs8);

    powers.push_back(powers1);
    powers.push_back(powers1);
    powers.push_back(powers1);
    powers.push_back(powers1);
    powers.push_back(powers2);
    powers.push_back(powers2);
    powers.push_back(powers2);
    powers.push_back(powers3);
    powers.push_back(powers3);
    powers.push_back(powers3);
    powers.push_back(powers4);
    powers.push_back(powers4);
    powers.push_back(powers4);
    powers.push_back(powers5);
    powers.push_back(powers6);
    powers.push_back(powers7);
    powers.push_back(powers8);
    powers.push_back(powers9);
    powers.push_back(powers10);

    angMom = 2;

    normalizeCoeffs();
}
