#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <armadillo>
#include "func.h"

using namespace std;
using namespace arma;


class Minimizer
{
public:
    Minimizer(Func* func);
    void solve();
    rowvec getMinPoint();
    double getMinValue();

private:
    mat m_X;
    rowvec m_fX;
    rowvec m_x0;
    Func* m_func;

    double m_alpha;
    double m_gamma;
    double m_rho;
    double m_sigma;
    double m_toler;
    double m_dim;

    void initializeSimplex();
    void advance();
    void sort();
    void setCentroid();
};

#endif // MINIMIZER_H
