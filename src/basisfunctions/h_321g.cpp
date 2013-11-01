#include "h_321g.h"

H_321G::H_321G()
{
    rowvec exp1 = {5.4471780, 0.8245470};
    rowvec exp2 = {0.1831920};

    rowvec coeffs1 = {0.1562850, 0.9046910};
    rowvec coeffs2 = {1.0};

    irowvec powers1 = {0, 0, 0};
    irowvec powers2 = {0, 0, 0};


    exponents.push_back(exp1);
    exponents.push_back(exp2);

    coeffs.push_back(coeffs1);
    coeffs.push_back(coeffs2);

    powers.push_back(powers1);
    powers.push_back(powers2);
}
