include(../defaults.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

LIBS += -L../src -lmyapp \
         -lyaml-cpp

OTHER_FILES += \
    configs/test.cfg
