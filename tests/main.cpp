#include <iostream>
#include <armadillo>
#include <unittest++/UnitTest++.h>
#include <src/integrator.h>

using namespace std;
using namespace arma;


TEST(OverlapIntegrals){
    Integrator integrator;
    rowvec3 RA = {1.2, 2.3, 3.4};
    integrator.setPositionA(RA);
    rowvec3 RB = {-1.3, 1.4, -2.4};
    integrator.setPositionB(RB);
    integrator.setAlpha(0.2);
    integrator.setBeta(0.3);
    integrator.setMaxAngMom(2);

    integrator.setE();
    CHECK_CLOSE(integrator.overlap(0, 0, 0, 0, 0, 0), 0.119172363581, 0.00001);
    CHECK_CLOSE(integrator.overlap(0, 0, 0, 0, 0, 2), 0.760605693318, 0.00001);
    CHECK_CLOSE(integrator.overlap(0, 0, 0, 2, 0, 0), 0.134617101901, 0.00001);
    CHECK_CLOSE(integrator.overlap(0, 0, 0, 2, 0, 2), 0.859180191173, 0.00001);
    CHECK_CLOSE(integrator.overlap(0, 0, 1, 0, 0, 0), -0.0643530763337, 0.00001);
    CHECK_CLOSE(integrator.overlap(0, 0, 1, 0, 0, 2), -0.410727074392, 0.00001);
    CHECK_CLOSE(integrator.overlap(0, 0, 1, 2, 0, 0), 0.0131108667517, 0.00001);
    CHECK_CLOSE(integrator.overlap(0, 0, 1, 2, 0, 2), 0.0836787959561, 0.00001);
    CHECK_CLOSE(integrator.overlap(1, 0, 0, 0, 0, 0), -0.178758545371, 0.00001);
    CHECK_CLOSE(integrator.overlap(1, 0, 0, 0, 0, 2), -1.14090853998, 0.00001);
    CHECK_CLOSE(integrator.overlap(1, 0, 0, 2, 0, 0), -0.201925652851, 0.00001);
    CHECK_CLOSE(integrator.overlap(1, 0, 0, 2, 0, 2), -1.28877028676, 0.00001);
    CHECK_CLOSE(integrator.overlap(1, 0, 1, 0, 0, 0), 0.0965296145005, 0.00001);
    CHECK_CLOSE(integrator.overlap(1, 0, 1, 0, 0, 2), 0.616090611588, 0.00001);
    CHECK_CLOSE(integrator.overlap(1, 0, 1, 2, 0, 0), -0.0196663001276, 0.00001);
    CHECK_CLOSE(integrator.overlap(1, 0, 1, 2, 0, 2), -0.125518193934, 0.00001);
}

int main()
{
    return UnitTest::RunAllTests();
}

