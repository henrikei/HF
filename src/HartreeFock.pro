TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo

QMAKE_CXXFLAGS += -std=c++0x
QMAKE_CXXFLAGS_RELEASE += -std=c++0x

SOURCES += main.cpp \
    hartreefock.cpp \
    system.cpp \
    primitivebasis.cpp \
    contractedbasis.cpp \
    integrator.cpp

HEADERS += \
    hartreefock.h \
    system.h \
    primitivebasis.h \
    contractedbasis.h \
    integrator.h

include(defaults.pri)
