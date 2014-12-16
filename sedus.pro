#/***************************************************************************
#**                                                                        **
#**  SEDUS, Segmental Duplication Simulator                                **
#**  Copyright (C) 2014 Diego A. Hartasánchez,Oriol Vallès-Codina,         **
#**  Marina Brasó-Vives, Juan Manuel Fuentes and Arcadi Navarro,           **
#**  Institut de Biologia Evolutiva UPF-CSIC                               **
#**                                                                        **
#**  This file is part of SEDUS.                                           **
#**                                                                        **
#**  SEDUS is free software: you can redistribute it and/or modify         **
#**  it under the terms of the GNU General Public License as published by  **
#**  the Free Software Foundation, either version 3 of the License, or     **
#**  (at your option) any later version.                                   **
#**                                                                        **
#**  SEDUS is distributed in the hope that it will be useful,              **
#**  but WITHOUT ANY WARRANTY; without even the implied warranty of        **
#**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         **
#**  GNU General Public License for more details.                          **
#**                                                                        **
#**  You should have received a copy of the GNU General Public License     **
#**  along with this program.  If not, see http://www.gnu.org/licenses/.   **
#**                                                                        **
#****************************************************************************
#**           Authors: Diego A. Hartasánchez,Oriol Vallès-Codina,          **
#**                   Marina Brasó-Vives, Juan Manuel Fuentes              **
#**                   and Arcadi Navarro                                   **
#**  Website/Contact: http://www.biologiaevolutiva.org/sedus/              **
#**             Date: 01.10.14                                             **
#**          Version: 1.10                                                 **
#****************************************************************************/

QT       += core gui



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
    pi_02_toFile_Sample.R \
    img/gnome-fs-directory-accept.png \
    img/header.png \
    img/logo_sedus.png \
    img/splash.png \
    img/start_64.png \
    img/stop.png \
    img/sedus.ico \
    COPYING \
    img/blocks.png \
    img/default.png

QMAKE_CXXFLAGS += -std=c++11

#copydata.commands = $(COPY_DIR) $$PWD/img $$OUT_PWD
#first.depends = $(first) copydata
#export(first.depends)
#export(copydata.commands)
#QMAKE_EXTRA_TARGETS += first copydata

RESOURCES += \
    sedus.qrc

win32:RC_ICONS += img/sedus.ico
