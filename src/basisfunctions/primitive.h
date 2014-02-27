#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <armadillo>

using namespace std;
using namespace arma;


class Primitive
{
public:
    Primitive(double exp, double coeff, irowvec pow, rowvec3 position);
    double getExp();
    double getCoeff();
    irowvec3 getPow();
    rowvec3 getPos();
private:
    double m_exp;
    double m_coeff;
    irowvec3 m_pow;
    rowvec3 m_position;
};

#endif // PRIMITIVE_H
