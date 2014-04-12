include(../defaults.pri)
CONFIG -= qt

TARGET = myapp
TEMPLATE = lib

SOURCES += \
    hartreefock/hartreefock.cpp \
    hartreefock/rhf.cpp \
    hartreefock/uhf.cpp \
    system/system.cpp \
    boysfunction/boysfunction.cpp \
    integrator/integrator.cpp \
    basisfunctions/primitive.cpp \
    basisfunctions/contracted.cpp \
    minimizer/minimizer.cpp \
    minimizer/func.cpp \
    minimizer/twodimtest.cpp \
    minimizer/hartreefockfunc.cpp \
    basisfunctions/basisfunctions.cpp \
    perturbation/mollerplesset.cpp \
    perturbation/rmp.cpp \
    perturbation/ump.cpp

HEADERS += \
    hartreefock/hartreefock.h \
    hartreefock/rhf.h \
    hartreefock/uhf.h \
    system/system.h \
    integrator/integrator.h \
    boysfunction/boysfunction.h \
    basisfunctions/primitive.h \
    basisfunctions/contracted.h \
    minimizer/minimizer.h \
    minimizer/func.h \
    minimizer/twodimtest.h \
    minimizer/hartreefockfunc.h \
    basisfunctions/basisfunctions.h \
    perturbation/mollerplesset.h \
    perturbation/rmp.h \
    perturbation/ump.h

OTHER_FILES += \
    ../inFiles/configFiles/H2O_431G.cfg \
    ../inFiles/configFiles/CH4_631Gs_UHF.cfg \
    ../inFiles/configFiles/CH4_631Gs.cfg \
    ../inFiles/configFiles/CH4_431G.cfg \
    ../inFiles/boys_tabulated.dat \
    ../inFiles/configFiles/CH4_631Gs_RHFvsUHF.cfg \
    ../inFiles/configFiles/H2_631Gss.cfg
