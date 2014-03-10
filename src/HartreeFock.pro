TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo -lboost_regex -lconfig++

QMAKE_CXXFLAGS += -std=c++0x
QMAKE_CXXFLAGS_RELEASE += -std=c++0x

SOURCES += main.cpp \
    hartreefock.cpp \
    system.cpp \
    boysfunction.cpp \
    basisfunctions/basisfunctions.cpp \
    integrator.cpp \
    basishandler.cpp \
    rhf.cpp \
    uhf.cpp \
    basisfunctions/primitive.cpp \
    basisfunctions/basisfunctions2.cpp \
    basisfunctions/contracted.cpp \
    minimizer/minimizer.cpp \
    minimizer/func.cpp \
    minimizer/twodimtest.cpp \
    minimizer/hartreefockfunc.cpp

HEADERS += \
    hartreefock.h \
    system.h \
    integrator.h \
    boysfunction.h \
    basisfunctions/basisfunctions.h \
    basishandler.h \
    rhf.h \
    uhf.h \
    basisfunctions/primitive.h \
    basisfunctions/basisfunctions2.h \
    basisfunctions/contracted.h \
    minimizer/minimizer.h \
    minimizer/func.h \
    minimizer/twodimtest.h \
    minimizer/hartreefockfunc.h

include(defaults.pri)

OTHER_FILES += \
    ../inFiles/configFiles/H2O_431G.cfg \
    ../inFiles/configFiles/CH4_631Gs_UHF.cfg \
    ../inFiles/configFiles/CH4_631Gs.cfg \
    ../inFiles/configFiles/CH4_431G.cfg \
    ../inFiles/boys_tabulated.dat \
    ../inFiles/configFiles/CH4_631Gs_RHFvsUHF.cfg \
    ../inFiles/configFiles/H2_631Gss.cfg
