#-------------------------------------------------
#
# Project created by QtCreator 2016-12-09T14:08:48
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Izhikevich_Model_Qt

TEMPLATE = app


SOURCES += main.cpp\
    qcustomplot.cpp \
    izhikevich_model_qt.cpp

HEADERS  += \
    qcustomplot.h \
    Neuron.h \
    izhikevich_model_qt.h

FORMS    += \
    izhikevich_model_qt.ui
