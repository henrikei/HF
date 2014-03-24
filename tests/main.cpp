#include <iostream>
#include <iomanip>
#include <armadillo>
#include <unittest++/UnitTest++.h>
#include <src/integrator.h>
#include <src/boysfunction.h>
#include <src/basisfunctions/basisfunctions2.h>
#include <src/basisfunctions/primitive.h>
#include <src/hartreefock.h>
#include <src/rhf.h>
#include <src/uhf.h>
#include <src/minimizer/minimizer.h>
#include <src/minimizer/func.h>
#include <src/minimizer/hartreefockfunc.h>
#include <src/minimizer/twodimtest.h>

using namespace std;
using namespace arma;


TEST(OverlapIntegrals){
    double alpha = 0.2;
    double beta = 0.3;
    double coeffA = 1.0;
    double coeffB = 1.0;
    irowvec3 powA = {0,0,0};
    irowvec3 powB = {0,0,0};
    rowvec3 RA = {1.2, 2.3, 3.4};
    rowvec3 RB = {-1.3, 1.4, -2.4};
    mat pos = zeros(2,3);
    pos.row(0) = RA;
    pos.row(1) = RB;
    Primitive *primitiveA = new Primitive(alpha, coeffA, powA, 0);
    Primitive *primitiveB = new Primitive(beta, coeffB, powB, 1);
    primitiveA->setPosPointer(&pos);
    primitiveB->setPosPointer(&pos);
    Integrator integrator(2);
    integrator.setPrimitiveA(primitiveA);
    integrator.setPrimitiveB(primitiveB);

    integrator.setE_AB("oneParticle");
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

    delete primitiveA;
    delete primitiveB;
}

TEST(KineticIntegrals){
    double alpha = 0.2;
    double beta = 0.3;
    double coeffA = 1.0;
    double coeffB = 1.0;
    irowvec3 powA = {1,1,1}; // Must be set equal or higher than the highest values below
    irowvec3 powB = {1,1,1}; // Must be set equal or higher than the highest values below
    rowvec3 RA = {1.2, 2.3, 3.4};
    rowvec3 RB = {-1.3, 1.4, -2.4};
    mat pos = zeros(2,3);
    pos.row(0) = RA;
    pos.row(1) = RB;
    Primitive *primitiveA = new Primitive(alpha, coeffA, powA, 0);
    Primitive *primitiveB = new Primitive(beta, coeffB, powB, 1);
    primitiveA->setPosPointer(&pos);
    primitiveB->setPosPointer(&pos);
    Integrator integrator(2);
    integrator.setPrimitiveA(primitiveA);
    integrator.setPrimitiveB(primitiveB);

    integrator.setE_AB("oneParticle");
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

    delete primitiveA;
    delete primitiveB;
}

TEST(BoysIntegrals){
    BoysFunction boys(3);

    boys.setx(2.3252);

    CHECK_CLOSE(boys.returnValue(1), 1.0007267355e-01, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(3), 2.5784878802e-02, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(6), 1.0688807154e-02, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(10), 5.8076406817e-03, 1.0E-10);

    boys.setx(8.9631);

    CHECK_CLOSE(boys.returnValue(1), 1.6505538584e-02, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(3), 7.6131457739e-04, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(6), 7.7859828141e-05, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(10), 1.9587800991e-05, 1.0E-10);

    boys.setx(17.6721);

    CHECK_CLOSE(boys.returnValue(1), 5.9646183419e-03, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(3), 7.1619859825e-05, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(6), 1.1232871899e-06, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(10), 4.4197332912e-08, 1.0E-10);

    boys.setx(31.1302);

    CHECK_CLOSE(boys.returnValue(1), 2.5511857399e-03, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(3), 9.8720995169e-06, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(6), 2.8346995237e-08, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(10), 1.1882122875e-10, 1.0E-10);

    boys.setx(44.8960);

    CHECK_CLOSE(boys.returnValue(1), 1.4730024561e-03, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(3), 2.7404345865e-06, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(6), 2.6232503507e-09, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(10), 2.5417117298e-12, 1.0E-10);

    boys.setx(57.1625);

    CHECK_CLOSE(boys.returnValue(1), 1.0252933095e-03, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(3), 1.1766761813e-06, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(6), 5.4571584880e-10, 1.0E-10);
    CHECK_CLOSE(boys.returnValue(10), 2.0120504238e-13, 1.0E-10);
}


TEST(Coulomb1){
    double alpha = 0.2;
    double beta = 0.3;
    double coeffA = 1.0;
    double coeffB = 1.0;
    irowvec3 powA = {0,0,0};
    irowvec3 powB = {0,0,0};
    rowvec3 RA = {1.2, 2.3, 3.4};
    rowvec3 RB = {-1.3, 1.4, -2.4};
    rowvec3 RC = {2.3, 0.9, 3.2};
    mat pos = zeros(2,3);
    pos.row(0) = RA;
    pos.row(1) = RB;
    Primitive *primitiveA = new Primitive(alpha, coeffA, powA, 0);
    Primitive *primitiveB = new Primitive(beta, coeffB, powB, 1);
    primitiveA->setPosPointer(&pos);
    primitiveB->setPosPointer(&pos);
    Integrator integrator(2);
    integrator.setPrimitiveA(primitiveA);
    integrator.setPrimitiveB(primitiveB);
    integrator.setNucleusPosition(RC);

    integrator.setE_AB("oneParticle");

    CHECK_CLOSE(2.788948987251e-02, integrator.coulomb1(0,0,0,0,0,0), 1e-5);
    CHECK_CLOSE(6.971203468743e-02, integrator.coulomb1(0,0,0,0,0,1), 1e-5);
    CHECK_CLOSE(2.024071525839e-01, integrator.coulomb1(0,0,0,0,0,2), 1e-5);
    CHECK_CLOSE(8.727033700014e-03, integrator.coulomb1(0,0,0,1,0,0), 1e-5);
    CHECK_CLOSE(2.134361291529e-02, integrator.coulomb1(0,0,0,1,0,1), 1e-5);
    CHECK_CLOSE(2.921666495443e-02, integrator.coulomb1(0,0,0,2,0,0), 1e-5);
    CHECK_CLOSE(3.185957751329e-02, integrator.coulomb1(0,1,0,0,0,0), 1e-5);
    CHECK_CLOSE(8.105746642202e-02, integrator.coulomb1(0,1,0,0,0,1), 1e-5);
    CHECK_CLOSE(9.596523510045e-03, integrator.coulomb1(0,1,0,1,0,0), 1e-5);
    CHECK_CLOSE(6.388444040338e-02, integrator.coulomb1(0,2,0,0,0,0), 1e-5);
    CHECK_CLOSE(-9.204700718547e-02, integrator.coulomb1(0,0,0,0,1,0), 1e-5);
    CHECK_CLOSE(-2.019226499202e-01, integrator.coulomb1(0,0,0,0,1,1), 1e-5);
    CHECK_CLOSE(-5.276399683274e-01, integrator.coulomb1(0,0,0,0,1,2), 1e-5);
    CHECK_CLOSE(-2.927318225377e-02, integrator.coulomb1(0,0,0,1,1,0), 1e-5);
    CHECK_CLOSE(-6.299787002318e-02, integrator.coulomb1(0,0,0,1,1,1), 1e-5);
    CHECK_CLOSE(-9.718105595370e-02, integrator.coulomb1(0,0,0,2,1,0), 1e-5);
    CHECK_CLOSE(-1.037280861539e-01, integrator.coulomb1(0,1,0,0,1,0), 1e-5);
    CHECK_CLOSE(-2.312309453843e-01, integrator.coulomb1(0,1,0,0,1,1), 1e-5);
    CHECK_CLOSE(-3.202910466576e-02, integrator.coulomb1(0,1,0,1,1,0), 1e-5);
    CHECK_CLOSE(-2.073449397904e-01, integrator.coulomb1(0,2,0,0,1,0), 1e-5);
    CHECK_CLOSE(3.319499900436e-01, integrator.coulomb1(0,0,0,0,2,0), 1e-5);
    CHECK_CLOSE(6.435114042344e-01, integrator.coulomb1(0,0,0,0,2,1), 1e-5);
    CHECK_CLOSE(1.536931448007e+00, integrator.coulomb1(0,0,0,0,2,2), 1e-5);
    CHECK_CLOSE(1.067865861209e-01, integrator.coulomb1(0,0,0,1,2,0), 1e-5);
    CHECK_CLOSE(2.033153544029e-01, integrator.coulomb1(0,0,0,1,2,1), 1e-5);
    CHECK_CLOSE(3.524622701603e-01, integrator.coulomb1(0,0,0,2,2,0), 1e-5);
//    CHECK_CLOSE(3.703919381580e-01, integrator.coulomb1(0,0,2,1,0,0), 1e-5);
//    CHECK_CLOSE(7.292169308884e-01, integrator.coulomb1(0,0,2,1,0,1), 1e-5);
//    CHECK_CLOSE(1.162963233448e-01, integrator.coulomb1(0,0,2,1,1,0), 1e-5);
//    CHECK_CLOSE(7.390872806284e-01, integrator.coulomb1(0,0,2,2,0,0), 1e-5);
//    CHECK_CLOSE(-1.637350724302e-02, integrator.coulomb1(0,1,0,0,0,0), 1e-5);
//    CHECK_CLOSE(-4.139721853567e-02, integrator.coulomb1(0,1,0,0,0,1), 1e-5);
//    CHECK_CLOSE(-1.213713540367e-01, integrator.coulomb1(0,1,0,0,0,2), 1e-5);
//    CHECK_CLOSE(2.136233458775e-02, integrator.coulomb1(0,1,0,0,1,0), 1e-5);
//    CHECK_CLOSE(5.306634838230e-02, integrator.coulomb1(0,1,0,0,1,1), 1e-5);
//    CHECK_CLOSE(-1.697963144320e-04, integrator.coulomb1(0,1,0,0,2,0), 1e-5);
//    CHECK_CLOSE(-1.907709578263e-02, integrator.coulomb1(0,1,0,1,0,0), 1e-5);
//    CHECK_CLOSE(-4.932098923684e-02, integrator.coulomb1(0,1,0,1,0,1), 1e-5);
//    CHECK_CLOSE(2.414126847830e-02, integrator.coulomb1(0,1,0,1,1,0), 1e-5);
//    CHECK_CLOSE(-3.842342257899e-02, integrator.coulomb1(0,1,0,2,0,0), 1e-5);
//    CHECK_CLOSE(5.356912442721e-02, integrator.coulomb1(0,1,1,0,0,0), 1e-5);
//    CHECK_CLOSE(1.187325132374e-01, integrator.coulomb1(0,1,1,0,0,1), 1e-5);
//    CHECK_CLOSE(3.128036711345e-01, integrator.coulomb1(0,1,1,0,0,2), 1e-5);
//    CHECK_CLOSE(-7.083519267815e-02, integrator.coulomb1(0,1,1,0,1,0), 1e-5);
//    CHECK_CLOSE(-1.544897767601e-01, integrator.coulomb1(0,1,1,0,1,1), 1e-5);
//    CHECK_CLOSE(-3.894393296797e-04, integrator.coulomb1(0,1,1,0,2,0), 1e-5);
//    CHECK_CLOSE(6.132616876012e-02, integrator.coulomb1(0,1,1,1,0,0), 1e-5);
//    CHECK_CLOSE(1.386353605834e-01, integrator.coulomb1(0,1,1,1,0,1), 1e-5);
//    CHECK_CLOSE(-7.911527548440e-02, integrator.coulomb1(0,1,1,1,1,0), 1e-5);
//    CHECK_CLOSE(1.229240947242e-01, integrator.coulomb1(0,1,1,2,0,0), 1e-5);
//    CHECK_CLOSE(3.609849112824e-02, integrator.coulomb1(0,2,0,0,0,0), 1e-5);
//    CHECK_CLOSE(9.032384498820e-02, integrator.coulomb1(0,2,0,0,0,1), 1e-5);
//    CHECK_CLOSE(2.625292648498e-01, integrator.coulomb1(0,2,0,0,0,2), 1e-5);
//    CHECK_CLOSE(-1.939589748931e-02, integrator.coulomb1(0,2,0,0,1,0), 1e-5);
//    CHECK_CLOSE(-4.913397190183e-02, integrator.coulomb1(0,2,0,0,1,1), 1e-5);
//    CHECK_CLOSE(6.878262296370e-02, integrator.coulomb1(0,2,0,0,2,0), 1e-5);
//    CHECK_CLOSE(4.131065513841e-02, integrator.coulomb1(0,2,0,1,0,0), 1e-5);
//    CHECK_CLOSE(1.052929737663e-01, integrator.coulomb1(0,2,0,1,0,1), 1e-5);
//    CHECK_CLOSE(-2.267402937768e-02, integrator.coulomb1(0,2,0,1,1,0), 1e-5);
//    CHECK_CLOSE(8.289710831960e-02, integrator.coulomb1(0,2,0,2,0,0), 1e-5);
//    CHECK_CLOSE(-3.786414780578e-02, integrator.coulomb1(1,0,0,0,0,0), 1e-5);
//    CHECK_CLOSE(-9.322262110550e-02, integrator.coulomb1(1,0,0,0,0,1), 1e-5);
//    CHECK_CLOSE(-2.671155215998e-01, integrator.coulomb1(1,0,0,0,0,2), 1e-5);
//    CHECK_CLOSE(-1.222106053447e-02, integrator.coulomb1(1,0,0,0,1,0), 1e-5);
//    CHECK_CLOSE(-2.972830178046e-02, integrator.coulomb1(1,0,0,0,1,1), 1e-5);
//    CHECK_CLOSE(-4.026352276293e-02, integrator.coulomb1(1,0,0,0,2,0), 1e-5);
//    CHECK_CLOSE(-1.576450345257e-02, integrator.coulomb1(1,0,0,1,0,0), 1e-5);
//    CHECK_CLOSE(-3.945885129414e-02, integrator.coulomb1(1,0,0,1,0,1), 1e-5);
//    CHECK_CLOSE(-4.918734877201e-03, integrator.coulomb1(1,0,0,1,1,0), 1e-5);
//    CHECK_CLOSE(-2.459437143524e-02, integrator.coulomb1(1,0,0,2,0,0), 1e-5);
//    CHECK_CLOSE(1.263894353489e-01, integrator.coulomb1(1,0,1,0,0,0), 1e-5);
//    CHECK_CLOSE(2.735756798558e-01, integrator.coulomb1(1,0,1,0,0,1), 1e-5);
//    CHECK_CLOSE(7.071773603054e-01, integrator.coulomb1(1,0,1,0,0,2), 1e-5);
//    CHECK_CLOSE(4.115384967585e-02, integrator.coulomb1(1,0,1,0,1,0), 1e-5);
//    CHECK_CLOSE(8.802219023191e-02, integrator.coulomb1(1,0,1,0,1,1), 1e-5);
//    CHECK_CLOSE(1.350111738398e-01, integrator.coulomb1(1,0,1,0,2,0), 1e-5);
//    CHECK_CLOSE(5.197526800952e-02, integrator.coulomb1(1,0,1,1,0,0), 1e-5);
//    CHECK_CLOSE(1.145639876363e-01, integrator.coulomb1(1,0,1,1,0,1), 1e-5);
//    CHECK_CLOSE(1.638641395940e-02, integrator.coulomb1(1,0,1,1,1,0), 1e-5);
//    CHECK_CLOSE(8.278875254192e-02, integrator.coulomb1(1,0,1,2,0,0), 1e-5);
//    CHECK_CLOSE(2.185667243163e-02, integrator.coulomb1(1,1,0,0,0,0), 1e-5);
//    CHECK_CLOSE(5.417205698627e-02, integrator.coulomb1(1,1,0,0,0,1), 1e-5);
//    CHECK_CLOSE(1.560020091608e-01, integrator.coulomb1(1,1,0,0,0,2), 1e-5);
//    CHECK_CLOSE(-2.926456829930e-02, integrator.coulomb1(1,1,0,0,1,0), 1e-5);
//    CHECK_CLOSE(-7.176178735649e-02, integrator.coulomb1(1,1,0,0,1,1), 1e-5);
//    CHECK_CLOSE(-5.223967979758e-04, integrator.coulomb1(1,1,0,0,2,0), 1e-5);
//    CHECK_CLOSE(9.269318129877e-03, integrator.coulomb1(1,1,0,1,0,0), 1e-5);
//    CHECK_CLOSE(2.337071697343e-02, integrator.coulomb1(1,1,0,1,0,1), 1e-5);
//    CHECK_CLOSE(-1.203714316117e-02, integrator.coulomb1(1,1,0,1,1,0), 1e-5);
//    CHECK_CLOSE(1.401501778682e-02, integrator.coulomb1(1,1,0,2,0,0), 1e-5);
//    CHECK_CLOSE(7.889586550718e-02, integrator.coulomb1(2,0,0,0,0,0), 1e-5);
//    CHECK_CLOSE(1.935977010010e-01, integrator.coulomb1(2,0,0,0,0,1), 1e-5);
//    CHECK_CLOSE(5.534914541236e-01, integrator.coulomb1(2,0,0,0,0,2), 1e-5);
//    CHECK_CLOSE(2.563391673303e-02, integrator.coulomb1(2,0,0,0,1,0), 1e-5);
//    CHECK_CLOSE(6.217850538435e-02, integrator.coulomb1(2,0,0,0,1,1), 1e-5);
//    CHECK_CLOSE(8.419480232293e-02, integrator.coulomb1(2,0,0,0,2,0), 1e-5);
//    CHECK_CLOSE(1.481688684288e-02, integrator.coulomb1(2,0,0,1,0,0), 1e-5);
//    CHECK_CLOSE(3.878852644576e-02, integrator.coulomb1(2,0,0,1,0,1), 1e-5);
//    CHECK_CLOSE(4.176920693786e-03, integrator.coulomb1(2,0,0,1,1,0), 1e-5);
//    CHECK_CLOSE(6.422210627967e-02, integrator.coulomb1(2,0,0,2,0,0), 1e-5);
    delete primitiveA;
    delete primitiveB;
}

TEST(H20_431G_RHF){
    // Simulation of H20 molecule with 4-31G basis.
    // Bond length (1.797) taken from Peter Atkins' book "Molecular Quantum Mechanics".
    // Bond angle set equal to 104.45 degrees.
    // Energy given by Atkins: -75.907.

    rowvec posO = {0.0, 0.0, 0.0};
    rowvec posH1 = {1.797, 0.0, 0.0};
    rowvec posH2 = {-0.448, 1.740, 0.0};
    rowvec charges = {8.0, 1.0, 1.0};
    int nElectrons = 10;

    mat nucleiPositions = zeros<mat>(3,3);
    nucleiPositions.row(0) = posO;
    nucleiPositions.row(1) = posH1;
    nucleiPositions.row(2) = posH2;

    BasisFunctions2 *basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/O_431G.dat", 0);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 1);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 2);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    RHF solver(system);
    solver.solve();
    double energy = solver.getEnergy();

    CHECK_CLOSE(-75.907, energy, 1.0E-3);

    delete basisFunctions;
    delete system;
}

TEST(H20_431G_UHF){
    // Simulation of H20 molecule with 4-31G basis.
    // Bond length (1.797) taken from Peter Atkins' book "Molecular Quantum Mechanics".
    // Bond angle set equal to 104.45 degrees.
    // Energy given by Atkins: -75.907.

    rowvec posO = {0.0, 0.0, 0.0};
    rowvec posH1 = {1.797, 0.0, 0.0};
    rowvec posH2 = {-0.448, 1.740, 0.0};
    rowvec charges = {8.0, 1.0, 1.0};
    int nElectrons = 10;

    mat nucleiPositions = zeros<mat>(3,3);
    nucleiPositions.row(0) = posO;
    nucleiPositions.row(1) = posH1;
    nucleiPositions.row(2) = posH2;

    BasisFunctions2 *basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/O_431G.dat", 0);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 1);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 2);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    UHF solver(system);
    solver.solve();
    double energy = solver.getEnergy();

    CHECK_CLOSE(-75.907, energy, 1.0E-3);

    delete basisFunctions;
    delete system;
}

TEST(CH4_431G_RHF){
    // Simulation of CH4 molecule.
    // Bond length (2.043) taken from Peter Atkins' book "Molecular Quantum Mechanics".
    // C-atom is placed at the origin. Direction vectors for the H-atoms are:
    // [1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]
    // Energy given by Atkins: -40.140.

    double d = 1.1795265999544056;
    rowvec posC = {0.0, 0.0, 0.0};
    rowvec posH1 = {d, d, d};
    rowvec posH2 = {-d, -d, d};
    rowvec posH3 = {d, -d, -d};
    rowvec posH4 = {-d, d, -d};
    rowvec charges = {6.0, 1.0, 1.0, 1.0, 1.0};
    int nElectrons = 10;

    mat nucleiPositions = zeros<mat>(5,3);
    nucleiPositions.row(0) = posC;
    nucleiPositions.row(1) = posH1;
    nucleiPositions.row(2) = posH2;
    nucleiPositions.row(3) = posH3;
    nucleiPositions.row(4) = posH4;

    BasisFunctions2 *basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/C_431G.dat", 0);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 1);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 2);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 3);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 4);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    RHF solver(system);
    solver.solve();
    double energy = solver.getEnergy();

    CHECK_CLOSE(energy, -40.140, 1.0E-3);

    delete basisFunctions;
    delete system;
}

TEST(CH4_431G_UHF){
    // Simulation of CH4 molecule.
    // Bond length (2.043) taken from Peter Atkins' book "Molecular Quantum Mechanics".
    // C-atom is placed at the origin. Direction vectors for the H-atoms are:
    // [1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]
    // Energy given by Atkins: -40.140.

    double d = 1.1795265999544056;
    rowvec posC = {0.0, 0.0, 0.0};
    rowvec posH1 = {d, d, d};
    rowvec posH2 = {-d, -d, d};
    rowvec posH3 = {d, -d, -d};
    rowvec posH4 = {-d, d, -d};
    rowvec charges = {6.0, 1.0, 1.0, 1.0, 1.0};
    int nElectrons = 10;

    mat nucleiPositions = zeros<mat>(5,3);
    nucleiPositions.row(0) = posC;
    nucleiPositions.row(1) = posH1;
    nucleiPositions.row(2) = posH2;
    nucleiPositions.row(3) = posH3;
    nucleiPositions.row(4) = posH4;

    BasisFunctions2 *basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/C_431G.dat", 0);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 1);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 2);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 3);
    basisFunctions->addContracteds("../inFiles/basisSets/H_431G.dat", 4);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    UHF solver(system);
    solver.solve();
    double energy = solver.getEnergy();

    CHECK_CLOSE(energy, -40.140, 1.0E-3);

    delete basisFunctions;
    delete system;
}

TEST(CH4_631Gss_RHF){
    // Simulation of CH4 molecule.
    // Bond length (2.043) taken from Peter Atkins' book "Molecular Quantum Mechanics".
    // C-atom is placed at the origin. Direction vectors for the H-atoms are:
    // [1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]
    // Energy given by Atkins: -40.202.

    double d = 1.1795265999544056;
    rowvec posC = {0.0, 0.0, 0.0};
    rowvec posH1 = {d, d, d};
    rowvec posH2 = {-d, -d, d};
    rowvec posH3 = {d, -d, -d};
    rowvec posH4 = {-d, d, -d};
    rowvec charges = {6.0, 1.0, 1.0, 1.0, 1.0};
    int nElectrons = 10;

    mat nucleiPositions = zeros<mat>(5,3);
    nucleiPositions.row(0) = posC;
    nucleiPositions.row(1) = posH1;
    nucleiPositions.row(2) = posH2;
    nucleiPositions.row(3) = posH3;
    nucleiPositions.row(4) = posH4;

    BasisFunctions2 *basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/C_631Gs.dat", 0);
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 1);
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 2);
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 3);
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 4);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    RHF solver(system);
    solver.solve();
    double energy = solver.getEnergy();

    CHECK_CLOSE(energy, -40.202, 1.0E-3);

    delete basisFunctions;
    delete system;
}

TEST(CH4_631Gss_UHF){
    // Simulation of CH4 molecule.
    // Bond length (2.043) taken from Peter Atkins' book "Molecular Quantum Mechanics".
    // C-atom is placed at the origin. Direction vectors for the H-atoms are:
    // [1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]
    // Energy given by Atkins: -40.202.

    double d = 1.1795265999544056;
    rowvec posC = {0.0, 0.0, 0.0};
    rowvec posH1 = {d, d, d};
    rowvec posH2 = {-d, -d, d};
    rowvec posH3 = {d, -d, -d};
    rowvec posH4 = {-d, d, -d};
    rowvec charges = {6.0, 1.0, 1.0, 1.0, 1.0};
    int nElectrons = 10;

    mat nucleiPositions = zeros<mat>(5,3);
    nucleiPositions.row(0) = posC;
    nucleiPositions.row(1) = posH1;
    nucleiPositions.row(2) = posH2;
    nucleiPositions.row(3) = posH3;
    nucleiPositions.row(4) = posH4;

    BasisFunctions2 *basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/C_631Gs.dat", 0);
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 1);
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 2);
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 3);
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 4);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    UHF solver(system);
    solver.solve();
    double energy = solver.getEnergy();

    CHECK_CLOSE(energy, -40.202, 1.0E-3);

    delete basisFunctions;
    delete system;
}

TEST(HF_631Gss_RHF){
    // Simulation of HF molecule.
    // Results checked against the following article:
    // "Full configuration interaction potential energy curves for breaking bonds
    // to hydrogen: An assessment of single-reference correlation methods"
    // by Dutta and Sherrill
    // JOURNAL OF CHEMICAL PHYSICS, VOLUME 118, NUMBER 4, 22 JANUARY 2003
    double d = 3.779451977;
    rowvec posH = {-d/2, 0.0, 0.0};
    rowvec posF = {d/2, 0.0, 0.0};
    rowvec charges = {1.0, 9.0};
    int nElectrons = 10;

    mat nucleiPositions = zeros<mat>(2,3);
    nucleiPositions.row(0) = posH;
    nucleiPositions.row(1) = posF;

    BasisFunctions2 *basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 0);
    basisFunctions->addContracteds("../inFiles/basisSets/F_631Gs.dat", 1);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    RHF solver(system);
    solver.solve();
    double energy = solver.getEnergy();

    CHECK_CLOSE(-99.747454, energy, 1.0E-6);

    delete basisFunctions;
    delete system;
}

TEST(HF_631Gss_UHF){
    // Simulation of HF molecule.
    // Results checked against the following article:
    // "Full configuration interaction potential energy curves for breaking bonds
    // to hydrogen: An assessment of single-reference correlation methods"
    // by Dutta and Sherrill
    // JOURNAL OF CHEMICAL PHYSICS, VOLUME 118, NUMBER 4, 22 JANUARY 2003
    double d = 3.779451977;
    rowvec posH = {-d/2, 0.0, 0.0};
    rowvec posF = {d/2, 0.0, 0.0};
    rowvec charges = {1.0, 9.0};
    int nElectrons = 10;

    mat nucleiPositions = zeros<mat>(2,3);
    nucleiPositions.row(0) = posH;
    nucleiPositions.row(1) = posF;

    BasisFunctions2 *basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/H_631Gss.dat", 0);
    basisFunctions->addContracteds("../inFiles/basisSets/F_631Gs.dat", 1);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    UHF solver(system);
    solver.solve();
    double energy = solver.getEnergy();

    CHECK_CLOSE(-99.865798, energy, 1.0E-6);

    delete basisFunctions;
    delete system;
}

TEST(FCl_RMP2_6311Gs){
    double d = 3.154519593;
    rowvec posF = {-0.5*d, 0.0, 0.0};
    rowvec posCl= {0.5*d, 0.0, 0.0};
    rowvec charges = {9.0, 17.0};
    int nElectrons = 26;

    mat nucleiPositions = zeros<mat>(2,3);
    nucleiPositions.row(0) = posF;
    nucleiPositions.row(1) = posCl;

    BasisFunctions2* basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/F_6311Gs.dat", 0);
    basisFunctions->addContracteds("../inFiles/basisSets/Cl_6311Gs.dat", 1);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    RHF solver(system,2);
    solver.solve();
    double energy = solver.getEnergy();

    CHECK_CLOSE(-559.293198, energy, 1.0E-6);

    delete basisFunctions;
    delete system;
}

TEST(FCl_UMP2_6311Gs){
    double d = 3.154519593;
    rowvec posF = {-0.5*d, 0.0, 0.0};
    rowvec posCl= {0.5*d, 0.0, 0.0};
    rowvec charges = {9.0, 17.0};
    int nElectrons = 26;

    mat nucleiPositions = zeros<mat>(2,3);
    nucleiPositions.row(0) = posF;
    nucleiPositions.row(1) = posCl;

    BasisFunctions2* basisFunctions = new BasisFunctions2;
    basisFunctions->addContracteds("../inFiles/basisSets/F_6311Gs.dat", 0);
    basisFunctions->addContracteds("../inFiles/basisSets/Cl_6311Gs.dat", 1);

    System *system;
    system = new System(basisFunctions, nucleiPositions, charges, nElectrons);

    UHF solver(system,2);
    solver.solve();
    double energy = solver.getEnergy();

    CHECK_CLOSE(-559.293198, energy, 1.0E-6);

    delete basisFunctions;
    delete system;
}

TEST(MINIMIZER_RBANANA){
    rowvec x = {-1.2, 1.0};
    Func *func = new TwoDimTest(x);
    Minimizer *minimizer= new Minimizer(func);
    minimizer->solve();
    x = minimizer->getMinPoint();

    CHECK_CLOSE(2.4, x(0), 1.0E-4);
    CHECK_CLOSE(5.76, x(1), 1.0E-4);

    delete func;
    delete minimizer;
}

int main()
{
    return UnitTest::RunAllTests();
}


