#include "density.h"

Density::Density(BasisFunctions *basisFunctions, field<mat> P, rowvec3 R1, rowvec3 R2,
                 double dx, double dy, double dz)
{
    m_basisFunctions = basisFunctions;
    m_R1 = R1, m_R2 = R2;
    m_dx = dx; m_dy = dy; m_dz = dz;
    m_P = P;
}

void Density::printDensity(string filename)
{
    if (m_P.n_elem == 1){ // Restricted determinant
        printDensity(m_P(0), filename);
    } else if (m_P.n_elem == 2){ // Unrestricted determinant
        mat Ptot = m_P(0) + m_P(1);
        printDensity(Ptot, filename);
    } else {
        cout << "Error in class Density: Dimension of field<mat> m_density invalid." << endl;
        exit(EXIT_FAILURE);
    }
}

void Density::printDensity(mat P, string filename)
{
    rowvec3 pos;
    int nBasisFunc = P.n_cols;
    Contracted *contractedA, *contractedB;

    int nx = (m_R2(0) - m_R1(0))/m_dx + 1;
    int ny = (m_R2(1) - m_R1(1))/m_dy + 1;
    int nz = (m_R2(2) - m_R1(2))/m_dz + 1;

    ofstream outfile;
    outfile.open(filename);
    outfile << m_R1(0) << endl; outfile << m_R1(1) << endl; outfile << m_R1(2) << endl;
    outfile << m_R2(0) << endl; outfile << m_R2(1) << endl; outfile << m_R2(2) << endl;
    outfile << m_dx << endl; outfile << m_dy << endl; outfile << m_dz << endl;
    outfile << nx << endl; outfile << ny << endl; outfile << nz << endl;

    double value = 0;
    for (int i = 0; i < nx; i++){
        pos(0) = m_R1(0) + i*m_dx;
        for (int j = 0; j < ny; j++){
            pos(1) = m_R1(1) + j*m_dy;
            for (int k = 0; k < nz; k++){
                pos(2) = m_R1(2) + k*m_dz;
                value = 0;

                for (int mu = 0; mu < nBasisFunc; mu++){
                    contractedA = m_basisFunctions->getContracted(mu);
                    for (int nu = 0; nu < nBasisFunc; nu++){
                        contractedB = m_basisFunctions->getContracted(nu);

                        value += P(mu,nu)*contractedA->getValue(pos)*contractedB->getValue(pos);
                    }
                }
                outfile << value << "  ";
            }
            outfile << endl;
        }
        outfile << endl;
    }
    outfile.close();
}
