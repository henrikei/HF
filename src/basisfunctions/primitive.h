#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <armadillo>

using namespace std;
using namespace arma;


class Primitive
{
public:
    Primitive(double exp, double coeff, irowvec pow, int posNum);
    double getExp();
    double getCoeff();
    irowvec3 getPow();
    rowvec3 getPos();
    void setPosPointer(mat* nucleiPositions);

private:
    double m_exp;
    double m_coeff;
    irowvec3 m_pow;
    int m_posNum;
    mat* m_nucleiPositions;
};

#endif // PRIMITIVE_H
