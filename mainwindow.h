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




#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QPushButton>
#include <QThread>
#include "sedus.h"
#include "qcustomplot.h"


//extern parameters params;

typedef QVector<QVector<double>> qvdouble;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:


    explicit MainWindow(QWidget *parent = 0);
    Ui::MainWindow *ui;
    ~MainWindow();

private slots:
    void handleButton();
    void stop();
    void RestDefault();
    void dirButton();
    void centerWidget();
    void setLog(const QString &value);
    void setBar(int value);
    void setChart(const qvdouble &x,const qvdouble &y);
    void setFixLinearText(int value);
    void setBurnin(int value);
    void setTotal(int value);
    void setFixLinearSliderVisible();
    void changeNdependecies(const QString &value);
//    void changeGenerationsDependecies(const QString &value);
    void changeR(const QString &value);
    void quitWork();
    void hotspots(int value);
    void setMinimumE1(double value);
    void setMaximumE1(double value);
    void setMinimumE2(double value);
    void setMaximumE2(double value);
    void setMinimumE3(double value);
    void setMaximumE3(double value);
    void setMinimumE4(double value);
    void setMaximumE4(double value);
    void setMinimumE5(double value);
    void setMaximumE5(double value);
    void changeMinimum(double value, int item);
    void changeMaximum(double value, int item);
    void changeAlpha1(double value);
    void changeAlpha2(double value);
    void changeAlpha3(double value);
    void changeAlpha4(double value);
    void changeAlpha5(double value);

signals:
    void setMinimum(double value, int item);
    void setMaximum(double value, int item);
private:
    QThread *thread;
    sedus *worker;
    bool checkRateHotSopts();
    bool checkGenerations();
    bool checkN();
    bool uidelete;
    double sumratiohots;
    int initmuN;
    //QPushButton *m_button;
};

#endif // MAINWINDOW_H
