include(../src/defaults.pri)

TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp

LIBS += -lunittest++

SOURCES += $$system(find $$SRC_DIR -name \'*.cpp\')
SOURCES = $$replace(SOURCES, $$SRC_DIR/main.cpp, )
message($$SOURCES)
