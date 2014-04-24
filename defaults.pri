LIBS += -larmadillo -llapack -lblas -lboost_regex -lconfig++\
        -larmadillo -lboost_mpi  -lboost_serialization -lboost_filesystem -lboost_system

#!nompi {
    # MPI Settings
    QMAKE_CXX = mpicxx
    QMAKE_CXX_RELEASE = $$QMAKE_CXX
    QMAKE_CXX_DEBUG = $$QMAKE_CXX
    QMAKE_LINK = $$QMAKE_CXX
    QMAKE_CC = mpicc

    QMAKE_CFLAGS += $$system(mpicc --showme:compile)
    QMAKE_LFLAGS += $$system(mpicxx --showme:link)
    QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
    QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

    DEFINES += USE_MPI
    #######
#}

COMMON_CXXFLAGS = -std=c++0x
QMAKE_CXXFLAGS += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG += $$COMMON_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$COMMON_CXXFLAGS

CURRENT_COMPILER = $$QMAKE_CXX
QMAKE_CXX = ccache $$CURRENT_COMPILER


INCLUDEPATH += $$PWD/src
SRC_DIR = $$PWD
