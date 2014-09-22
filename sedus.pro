#-------------------------------------------------
#
# Project created by QtCreator 2014-09-01T12:21:58
#
#-------------------------------------------------

QT       += core gui

QXT = core

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = sedus
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    sedus.cpp \
    qcustomplot.cpp

HEADERS  += mainwindow.h \
    sedus.h \
    qcustomplot.h

FORMS    += mainwindow.ui

OTHER_FILES += \
    pi_02_toFile_Sample.R

QMAKE_CXXFLAGS += -std=c++11

copydata.commands = $(COPY_DIR) $$PWD/img $$OUT_PWD
first.depends = $(first) copydata
export(first.depends)
export(copydata.commands)
QMAKE_EXTRA_TARGETS += first copydata
