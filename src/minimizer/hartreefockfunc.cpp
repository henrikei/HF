#include "hartreefockfunc.h"

HartreeFockFunc::HartreeFockFunc(MollerPlesset *solver, System *system)
{
    m_solver = solver;
    m_system = system;
    m_nucleiPositions = system->getNucleiPositions();
    m_nucleiPositionsTransformed.copy_size(m_nucleiPositions);
    m_nucleiPositionsTransformed.fill(0.0);

    int N = m_nucleiPositions.n_rows;

    m_trans = zeros(N, 3);
    m_rotx = zeros(3,3);
    m_roty = zeros(3,3);
    m_rotz = zeros(3,3);

    if (N < 2){
        cout << "Error: No degree of freedom to vary during minimization" << endl;
        exit(EXIT_FAILURE);
    } else if (N > 1){
        m_x = zeros<rowvec>(1);

        pair<int,int> a(1,0);
        m_map.push_back(a);
    }
    if (N > 2){
        m_x = zeros<rowvec>(3*N - 6);

        pair<int,int> a(2,0);
        pair<int,int> b(2,1);
        m_map.push_back(a); m_map.push_back(b);
        int row, col;
        for (uint i = 0; i < m_x.n_elem; i++){
            row = 3 + i/3;
            col = i % 3;
            pair<int,int> c(row,col);
            m_map.push_back(c);
        }
    }

    calcTransfMatrices();
    transfPositions();
    matToVec();
}


rowvec HartreeFockFunc::getx()
{
    return m_x;
}


double HartreeFockFunc::getValue(rowvec x)
{
    m_x = x;
    vecToMat();

    int my_rank = 0;
#ifdef RUN_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#endif
    if (my_rank == 0){
        cout << m_nucleiPositionsTransformed << endl;
    }

    m_system->setNucleiPositions(m_nucleiPositionsTransformed);
    m_solver->solve();
    return m_solver->getEnergy();
}


mat HartreeFockFunc::getNucleiPositions()
{
    vecToMat();
    transfPositionsInverse();
    return m_nucleiPositions;
}


// Calculates on the basis of m_nucleiPositions the matrices necessary to place:
// 1) atom 1 at the origin
// 2) atom 2 at the x-axis
// 3) atom 3 in the xy-plane.
void HartreeFockFunc::calcTransfMatrices()
{
    // Define translation which moves particle 1 to the origin
    for (uint i = 0; i < m_trans.n_rows; i++){
        m_trans.row(i) = -m_nucleiPositions.row(0);
    }
    mat nucleiPositions = m_nucleiPositions;
    nucleiPositions = nucleiPositions + m_trans;

    // Define necessary rotations to place atom 2 at the x-axis.
    // cx = cos(theta_x) etc
    double cz = nucleiPositions(1,0)/sqrt(nucleiPositions(1,0)*nucleiPositions(1,0) + nucleiPositions(1,1)*nucleiPositions(1,1));
    double sz = nucleiPositions(1,1)/sqrt(nucleiPositions(1,0)*nucleiPositions(1,0) + nucleiPositions(1,1)*nucleiPositions(1,1));

    m_rotz(0,0) = cz; m_rotz(0,1) = sz;
    m_rotz(1,0) = -sz; m_rotz(1,1) = cz;
    m_rotz(2,2) = 1.0;

    nucleiPositions = nucleiPositions*m_rotz.t();

    double cy = nucleiPositions(1,0)/sqrt(nucleiPositions(1,0)*nucleiPositions(1,0) + nucleiPositions(1,2)*nucleiPositions(1,2));
    double sy = nucleiPositions(1,2)/sqrt(nucleiPositions(1,0)*nucleiPositions(1,0) + nucleiPositions(1,2)*nucleiPositions(1,2));

    m_roty(0,0) = cy; m_roty(0,2) = sy;
    m_roty(1,1) = 1.0;
    m_roty(2,0) = -sy; m_roty(2,2) = cy;

    nucleiPositions = nucleiPositions*m_roty.t();

    // Define necessary rotations to place atom 3 in the xy-plane
    if (nucleiPositions.n_rows > 2){
        // If atom 3 is already in the xy-plane, the rotation matrix is the identity
        if (nucleiPositions(2,2) == 0){
            m_rotx(0,0) = 1.0; m_rotx(1,1) = 1.0; m_rotx(2,2) = 1.0;
        } else {
            double cx = nucleiPositions(2,1)/sqrt(nucleiPositions(2,1)*nucleiPositions(2,1) + nucleiPositions(2,2)*nucleiPositions(2,2));
            double sx = nucleiPositions(2,2)/sqrt(nucleiPositions(2,1)*nucleiPositions(2,1) + nucleiPositions(2,2)*nucleiPositions(2,2));

            m_rotx(0,0) = 1.0;
            m_rotx(1,1) = cx; m_rotx(1,2) = sx;
            m_rotx(2,1) = -sx; m_rotx(2,2) = cx;
        }
    }
}

void HartreeFockFunc::transfPositions()
{
    m_nucleiPositionsTransformed = m_nucleiPositions + m_trans;
    if (m_nucleiPositionsTransformed.n_rows > 2){
        m_nucleiPositionsTransformed = m_nucleiPositions*m_rotz.t()*m_roty.t()*m_rotx.t();
    } else {
        m_nucleiPositionsTransformed = m_nucleiPositions*m_rotz.t()*m_roty.t();
    }
}

void HartreeFockFunc::transfPositionsInverse()
{
    if (m_nucleiPositions.n_rows > 2){
        m_nucleiPositions = m_nucleiPositionsTransformed*inv(m_rotx.t())*inv(m_roty.t())*inv(m_rotz.t());
    } else {
        m_nucleiPositions = m_nucleiPositionsTransformed*inv(m_roty.t())*inv(m_rotz.t());
    }
    m_nucleiPositions = m_nucleiPositionsTransformed - m_trans;
}

void HartreeFockFunc::matToVec()
{
//    m_x(0) = m_nucleiPositionsTransformed(1,0);

//    if (m_nucleiPositionsTransformed.n_rows > 2){
//        m_x(1) = m_nucleiPositionsTransformed(2,0);
//        m_x(2) = m_nucleiPositionsTransformed(2,1);
//        int row, col;
//        for (uint i = 0; i < m_x.n_elem-3; i++){
//            row = 3 + i/3;
//            col = i % 3;
//            m_x(3+i) = m_nucleiPositionsTransformed(row,col);
//        }
//    }
    int row, col;
    for(uint i = 0; i < m_x.n_elem; i++){
        row = m_map.at(i).first;
        col = m_map.at(i).second;
        m_x(i) = m_nucleiPositionsTransformed(row,col);
    }
}

void HartreeFockFunc::vecToMat()
{
    int row, col;
    for (uint i = 0; i < m_x.n_elem; i++){
        row = m_map.at(i).first;
        col = m_map.at(i).second;
        m_nucleiPositionsTransformed(row,col) = m_x(i);
    }
//    m_nucleiPositionsTransformed(1,0) = m_x(0);
//    if (m_nucleiPositionsTransformed.n_rows > 2){
//        m_nucleiPositionsTransformed(2,0) = m_x(1);
//        m_nucleiPositionsTransformed(2,1) = m_x(2);
//        int row, col;
//        for (uint i = 0; i < m_x.n_elem-3; i++){
//            row = 3 + i/3;
//            col = i % 3;
//            m_nucleiPositionsTransformed(row,col) = m_x(3+i);
//        }
//    }
}
