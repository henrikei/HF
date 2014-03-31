include(../defaults.pri)
CONFIG -= qt

TARGET = myapp
TEMPLATE = lib

SOURCES += \
    hartreefock.cpp \
    system.cpp \
    boysfunction.cpp \
    integrator.cpp \
    rhf.cpp \
    uhf.cpp \
    basisfunctions/primitive.cpp \
    basisfunctions/basisfunctions2.cpp \
    basisfunctions/contracted.cpp \
    minimizer/minimizer.cpp \
    minimizer/func.cpp \
    minimizer/twodimtest.cpp \
    minimizer/hartreefockfunc.cpp \
    perturbation/mollerplessetpt.cpp \
    perturbation/restrictedmollerplessetpt.cpp \
    perturbation/unrestrictedmollerplessetpt.cpp

HEADERS += \
    hartreefock.h \
    system.h \
    integrator.h \
    boysfunction.h \
    rhf.h \
    uhf.h \
    basisfunctions/primitive.h \
    basisfunctions/basisfunctions2.h \
    basisfunctions/contracted.h \
    minimizer/minimizer.h \
    minimizer/func.h \
    minimizer/twodimtest.h \
    minimizer/hartreefockfunc.h \
    perturbation/mollerplessetpt.h \
    perturbation/restrictedmollerplessetpt.h \
    perturbation/unrestrictedmollerplessetpt.h

OTHER_FILES += \
    ../inFiles/configFiles/H2O_431G.cfg \
    ../inFiles/configFiles/CH4_631Gs_UHF.cfg \
    ../inFiles/configFiles/CH4_631Gs.cfg \
    ../inFiles/configFiles/CH4_431G.cfg \
    ../inFiles/boys_tabulated.dat \
    ../inFiles/configFiles/CH4_631Gs_RHFvsUHF.cfg \
    ../inFiles/configFiles/H2_631Gss.cfg
