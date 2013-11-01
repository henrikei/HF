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
    basisfunctions.cpp \
    basisfunctions/basisfunctions.cpp \
    basisfunctions/h_321g.cpp \
    integrator.cpp

HEADERS += \
    hartreefock.h \
    system.h \
    integrator.h \
    boysfunction.h \
    basisfunctions/basisfunctions.h \
    basisfunctions/h_321g.h

include(defaults.pri)
