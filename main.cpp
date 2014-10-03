/***************************************************************************
**                                                                        **
**  SEDUS, Segmental Duplication Simulator                                **
**  Copyright (C) 2014 Diego A. Hartasánchez, Oriol Vallès-Codina,        **
**  Marina Brasó-Vives, Juan Manuel Fuentes and Arcadi Navarro,           **
**  Institut de Biologia Evolutiva UPF-CSIC                               **
**                                                                        **
**  This file is part of SEDUS.                                           **
**                                                                        **
**  SEDUS is free software: you can redistribute it and/or modify         **
**  it under the terms of the GNU General Public License as published by  **
**  the Free Software Foundation, either version 3 of the License, or     **
**  (at your option) any later version.                                   **
**                                                                        **
**  SEDUS is distributed in the hope that it will be useful,              **
**  but WITHOUT ANY WARRANTY; without even the implied warranty of        **
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         **
**  GNU General Public License for more details.                          **
**                                                                        **
**  You should have received a copy of the GNU General Public License     **
**  along with this program.  If not, see http://www.gnu.org/licenses/.   **
**                                                                        **
****************************************************************************
**           Authors: Diego A. Hartasánchez, Oriol Vallès-Codina,         **
**                   Marina Brasó-Vives, Juan Manuel Fuentes              **
**                   and Arcadi Navarro                                   **
**  Website/Contact: http://www.biologiaevolutiva.org/sedus/              **
**             Date: 01.10.14                                             **
**          Version: 1.10                                                 **
****************************************************************************/

#include "mainwindow.h"
#include <QApplication>
#include <QSplashScreen>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    #ifndef __APPLE__
        a.setWindowIcon(QIcon(":/img/sedus.ico"));
    #endif
    QPixmap pixmap(":/img/splash.png");
    QSplashScreen splash(pixmap);
    splash.show();
    splash.showMessage(QObject::tr("Initiating SEDUS..."),
                        Qt::AlignLeft | Qt::AlignTop, Qt::black);
    a.processEvents();
    QThread::sleep(2);

    splash.showMessage(QObject::tr("GNU GENERAL PUBLIC LICENSE, http://www.gnu.org/licenses/gpl.txt"),
                        Qt::AlignLeft | Qt::AlignTop, Qt::black);

    a.processEvents();
    QThread::sleep(3);
    MainWindow w;
    w.show();
    splash.finish(&w);
    splash.raise();
    return a.exec();
}
