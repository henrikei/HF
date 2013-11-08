TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo

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
    basisfunctions/n_431g.cpp

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
    basisfunctions/n_431g.h

include(defaults.pri)
