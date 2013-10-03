TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo

SOURCES += main.cpp \
    hartreefock.cpp

HEADERS += \
    hartreefock.h

