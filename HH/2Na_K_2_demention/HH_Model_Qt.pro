#-------------------------------------------------
#
# Project created by QtCreator 2016-12-09T14:08:48
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = HH_Model_Qt
TEMPLATE = app


SOURCES += main.cpp\
        hh_model_qt.cpp \
    qcustomplot.cpp

HEADERS  += hh_model_qt.h \
    qcustomplot.h \
    Neuron.h

FORMS    += hh_model_qt.ui
