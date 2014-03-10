#include "minimizer.h"

Minimizer::Minimizer(Func* func, rowvec x)
{
    m_func = func;
    m_dim = x.n_elem;
    m_x0 = zeros<rowvec>(m_dim);
    initializeSimplex(x);

    m_alpha = 1.0;
    m_gamma = 2.0;
    m_rho = -0.5;
    m_sigma = 0.5;
    m_toler = 1.0E-8;
}

rowvec Minimizer::solve()
{
    while (fabs(m_fX(m_dim) - m_fX(0)) > m_toler){
        advance();
    }
    return m_X.row(0);
}

void Minimizer::initializeSimplex(rowvec x)
{
    m_X.set_size(m_dim + 1, m_dim);
    m_X.row(0) = x;

    double initialStep = 1.0;
    rowvec xtemp = zeros<rowvec>(m_dim);

    for (uint i = 1; i < m_dim+1; i++){
        xtemp = x;
        xtemp(i-1) += initialStep;
        m_X.row(i) = xtemp;
    }

    m_fX = zeros<rowvec>(m_dim+1);
    for (uint i = 0; i < m_dim+1; i++){
        m_fX(i) = m_func->getValue(m_X.row(i));
    }

    sort();
    setCentroid();
}

void Minimizer::advance()
{
    rowvec xNew = m_x0 + m_alpha*(m_x0 - m_X.row(m_dim));
    double fxNew = m_func->getValue(xNew);

    if (m_fX(0) <= fxNew && fxNew < m_fX(m_dim-1)){
        m_X.row(m_dim) = xNew;
        m_fX(m_dim) = fxNew;
    } else if (fxNew < m_fX(0)){
        m_fX(m_dim) = fxNew;
        m_X.row(m_dim) = xNew;
        xNew = m_x0 + m_gamma*(m_x0 - m_X.row(m_dim));
        fxNew = m_func->getValue(xNew);
        if (fxNew < m_fX(m_dim)){
            m_X.row(m_dim) = xNew;
            m_fX(m_dim) = fxNew;
        }
    } else {
        xNew = m_x0 + m_rho*(m_x0 - m_X.row(m_dim));
        fxNew = m_func->getValue(xNew);
        if (fxNew < m_fX(m_dim)){
            m_X.row(m_dim) = xNew;
            m_fX(m_dim) = fxNew;
        } else {
            for (int i = 1; i < m_dim + 1; i++){
                m_X.row(i) = m_X.row(0) + m_sigma*(m_X.row(i) - m_X.row(0));
                m_fX(i) = m_func->getValue(m_X.row(i));
            }
        }
    }
    sort();
    setCentroid();
}

void Minimizer::sort()
{
    double f1, f2;
    rowvec x1, x2;
    bool unsorted = true;
    while (unsorted){
        unsorted = false;
        for (uint i = 1; i < m_dim+1; i++){
            x1 = m_X.row(i-1);
            x2 = m_X.row(i);
            f1 = m_fX(i-1);
            f2 = m_fX(i);
            if (f1 > f2){
                m_fX(i-1) = f2;
                m_fX(i) = f1;
                m_X.row(i-1) = x2;
                m_X.row(i) = x1;
                unsorted = true;
            }
        }
    }
}

void Minimizer::setCentroid()
{
    m_x0 = zeros<rowvec>(m_dim);
    for (uint i = 0; i < m_dim; i++){
        m_x0 += m_X.row(i);
    }
    m_x0 /= m_dim;
}
