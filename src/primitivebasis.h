#ifndef PRIMITIVEBASIS_H
#define PRIMITIVEBASIS_H

class PrimitiveBasis
{
public:
    PrimitiveBasis(double alpha, double d, char type);
    double getAlpha();
    double getCoeff();
private:
    double m_alpha;
    double m_coeff;
    char m_type;
};

#endif // PRIMITIVEBASIS_H
