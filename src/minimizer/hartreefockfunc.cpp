#include "hartreefockfunc.h"

HartreeFockFunc::HartreeFockFunc(HartreeFock *solver, System *system)
{
    m_solver = solver;
    m_system = system;

    m_trans = zeros(m_system->getNucleiPositions().n_rows, 3);
    m_rotx = zeros(3,3);
    m_roty = zeros(3,3);
    m_rotz = zeros(3,3);

    m_nucleiPositions = m_system->getNucleiPositions();
    if (m_nucleiPositions.n_rows < 2){
        cout << "Error: No degree of freedom to vary during minimization" << endl;
        exit(EXIT_FAILURE);
    }
}

rowvec HartreeFockFunc::getxInitial()
{
    calcTransfMatrices();
    transfPositions();
    vec x = matToVec(m_nucleiPositions);
    return x;
}

double HartreeFockFunc::getValue(rowvec x)
{
    m_nucleiPositions = vecToMat(x);
    m_system->setNucleiPositions(m_nucleiPositions);
    return m_solver->getEnergy();
}

void HartreeFockFunc::calcTransfMatrices()
{
    // Define necessary translation to place particle 1 at the origin
    for (uint i = 0; i < m_trans.n_rows; i++){
        m_trans.row(i) = -m_nucleiPositions.row(0);
    }
    // Define necessary rotations to place atom 2 at the x-axis.
    // cx = cos(theta_x) etc
    double cz = m_nucleiPositions(1,0)/sqrt(m_nucleiPositions(1,0)*m_nucleiPositions(1,0) + m_nucleiPositions(1,1)*m_nucleiPositions(1,1));
    double sz = m_nucleiPositions(1,1)/sqrt(m_nucleiPositions(1,0)*m_nucleiPositions(1,0) + m_nucleiPositions(1,1)*m_nucleiPositions(1,1));
    double cy = m_nucleiPositions(1,0)/sqrt(m_nucleiPositions(1,0)*m_nucleiPositions(1,0) + m_nucleiPositions(1,2)*m_nucleiPositions(1,2));
    double sy = m_nucleiPositions(1,2)/sqrt(m_nucleiPositions(1,0)*m_nucleiPositions(1,0) + m_nucleiPositions(1,2)*m_nucleiPositions(1,2));

    m_rotz(0,0) = cz; m_rotz(0,1) = sz;
    m_rotz(1,0) = -sz; m_rotz(1,1) = cz;
    m_roty(0,0) = cy; m_roty(0,2) = -sz;
    m_roty(2,0) = sy; m_roty(2,2) = cy;

    mat nucleiPositions = m_nucleiPositions;
    nucleiPositions = m_nucleiPositions + m_trans;
    nucleiPositions = m_nucleiPositions*m_rotz.t();
    nucleiPositions = m_nucleiPositions*m_roty.t();

    if (nucleiPositions.n_rows > 2){
        double cx = nucleiPositions(2,1)/sqrt(nucleiPositions(2,1)*nucleiPositions(2,1) + nucleiPositions(2,2)*nucleiPositions(2,2));
        double sx = nucleiPositions(2,2)/sqrt(nucleiPositions(2,1)*nucleiPositions(2,1) + nucleiPositions(2,2)*nucleiPositions(2,2));

        m_rotx(1,1) = cx; m_rotx(1,2) = -sx;
        m_rotx(2,1) = sx; m_rotx(2,2) = cx;
    }
}

void HartreeFockFunc::transfPositions()
{
    m_nucleiPositions = m_nucleiPositions + m_trans;
    m_nucleiPositions = m_nucleiPositions*m_rotz.t()*m_roty.t()*m_rotx.t();
}

void HartreeFockFunc::transfPositionsInverse()
{
    m_nucleiPositions = m_nucleiPositions*inv(m_rotx.t())*inv(m_roty.t())*inv(m_rotz.t());
    m_nucleiPositions = m_nucleiPositions - m_trans;
}

rowvec HartreeFockFunc::matToVec(mat nucleiPositions)
{
    int N = nucleiPositions.n_rows;
    rowvec x;
    if (N == 2){
       x = zeros<rowvec>(1);
       x = nucleiPositions(1,0);
    } else if (N > 2){
        x = zeros<rowvec>(2*N - 6);
        x(1) = m_nucleiPositions(2,0);
        x(2) = m_nucleiPositions(2,1);
        for (uint i = 3; i < x.n_elem; i++){
            x(i) = m_nucleiPositions(i+3);
        }
    }
    return x;
}

mat HartreeFockFunc::vecToMat(rowvec x)
{
    int N = m_system->getNucleiPositions().n_rows;
    mat m_nucleiPositions = zeros(N, 3);
    m_nucleiPositions(1,0) = x(0);
    if (m_nucleiPositions.n_rows > 2){
        m_nucleiPositions(2,0) = x(1);
        m_nucleiPositions(2,1) = x(2);
        for (uint i = 3; i < x.n_elem; i++){
            m_nucleiPositions(i+3) = x(i);
        }
    }
    return m_nucleiPositions;
}
