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
    basisfunctions/h_321g.cpp \
    integrator.cpp \
    basishandler.cpp \
    basisfunctions/h_theijssen.cpp \
    basisfunctions/o_321g.cpp \
    basisfunctions/h_431g.cpp \
    basisfunctions/n_431g.cpp \
    basisfunctions/o_431g.cpp \
    rhf.cpp \
    uhf.cpp \
    basisfunctions/h_6311gss.cpp \
    basisfunctions/o_6311gss.cpp \
    basisfunctions/o_6311g.cpp

HEADERS += \
    hartreefock.h \
    system.h \
    integrator.h \
    boysfunction.h \
    basisfunctions/basisfunctions.h \
    basisfunctions/h_321g.h \
    basishandler.h \
    basisfunctions/h_theijssen.h \
    basisfunctions/o_321g.h \
    basisfunctions/h_431g.h \
    basisfunctions/n_431g.h \
    basisfunctions/o_431g.h \
    rhf.h \
    uhf.h \
    basisfunctions/h_6311gss.h \
    basisfunctions/o_6311gss.h \
    basisfunctions/o_6311g.h

include(defaults.pri)

OTHER_FILES += \
    ../inFiles/configFiles/H2O_431G.cfg \
    ../inFiles/configFiles/CH4_631Gs_UHF.cfg \
    ../inFiles/configFiles/CH4_631Gs.cfg \
    ../inFiles/configFiles/CH4_431G.cfg \
    ../inFiles/boys_tabulated.dat
