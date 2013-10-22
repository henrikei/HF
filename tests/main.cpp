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

TEST(KineticIntegrals){
    Integrator integrator;
    rowvec3 RA = {1.2, 2.3, 3.4};
    integrator.setPositionA(RA);
    rowvec3 RB = {-1.3, 1.4, -2.4};
    integrator.setPositionB(RB);
    integrator.setAlpha(0.2);
    integrator.setBeta(0.3);
    integrator.setMaxAngMom(2);

    integrator.setE();
    CHECK_CLOSE(integrator.kinetic(0, 0, 0, 0, 0, 0), -0.0967870268058, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 0, 0, 0, 1), -0.158190730148, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 0, 0, 1, 0), 0.237286095222, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 0, 0, 1, 1), 0.251402082662, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 0, 1, 0, 0), -0.0245468374367, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 0, 1, 0, 1), -0.0330608009181, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 0, 1, 1, 0), 0.0495912013772, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 0, 1, 1, 1), 0.0176714824377, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 1, 0, 0, 0), 0.0368202561551, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 1, 0, 0, 1), 0.0495912013772, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 1, 0, 1, 0), -0.0743868020658, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 1, 0, 1, 1), -0.0265072236566, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 1, 1, 0, 0), -0.0604904731258, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 1, 1, 0, 1), -0.086882171055, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 1, 1, 1, 0), 0.130323256583, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 0, 1, 1, 1, 1), 0.0788748150527, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 0, 0, 0, 0), -0.0681856595464, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 0, 0, 0, 1), -0.0918355581059, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 0, 0, 1, 0), 0.137753337159, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 0, 0, 1, 1), 0.0490874512159, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 0, 1, 0, 0), -0.0142503452233, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 0, 1, 0, 1), -0.00917293898306, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 0, 1, 1, 0), 0.0137594084746, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 0, 1, 1, 1), -0.0551617848829, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 1, 0, 0, 0), 0.021375517835, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 1, 0, 0, 1), 0.0137594084746, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 1, 0, 1, 0), -0.0206391127119, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 1, 0, 1, 1), 0.0827426773243, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 1, 1, 0, 0), -0.0374492116616, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 1, 1, 0, 1), -0.0334264444581, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 1, 1, 1, 0), 0.0501396666872, 0.00001);
    CHECK_CLOSE(integrator.kinetic(0, 1, 1, 1, 1, 1), -0.0841098520403, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 0, 0, 0, 0), 0.10227848932, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 0, 0, 0, 1), 0.137753337159, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 0, 0, 1, 0), -0.206630005738, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 0, 0, 1, 1), -0.0736311768238, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 0, 1, 0, 0), 0.021375517835, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 0, 1, 0, 1), 0.0137594084746, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 0, 1, 1, 0), -0.0206391127119, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 0, 1, 1, 1), 0.0827426773243, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 1, 0, 0, 0), -0.0320632767525, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 1, 0, 0, 1), -0.0206391127119, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 1, 0, 1, 0), 0.0309586690678, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 1, 0, 1, 1), -0.124114015986, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 1, 1, 0, 0), 0.0561738174925, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 1, 1, 0, 1), 0.0501396666872, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 1, 1, 1, 0), -0.0752095000308, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 0, 1, 1, 1, 1), 0.126164778061, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 0, 0, 0, 0), -0.0088092211159, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 0, 0, 0, 1), -0.0536149790098, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 0, 0, 1, 0), 0.0804224685147, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 0, 0, 1, 1), 0.278928221561, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 0, 1, 0, 0), -0.00831956570842, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 0, 1, 0, 1), -0.0312453234111, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 0, 1, 1, 0), 0.0468679851166, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 0, 1, 1, 1), 0.136830793422, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 1, 0, 0, 0), 0.0124793485626, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 1, 0, 0, 1), 0.0468679851166, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 1, 0, 1, 0), -0.0703019776749, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 1, 0, 1, 1), -0.205246190134, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 1, 1, 0, 0), -0.0158372863654, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 1, 1, 0, 1), -0.0634703676663, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 1, 1, 1, 0), 0.0952055514994, 0.00001);
    CHECK_CLOSE(integrator.kinetic(1, 1, 1, 1, 1, 1), 0.28653192666, 0.00001);
}

int main()
{
    return UnitTest::RunAllTests();
}

