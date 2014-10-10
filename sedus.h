/***************************************************************************
**                                                                        **
**  SEDUS, Segmental Duplication Simulator                                **
**  Copyright (C) 2014 Diego A. Hartasánchez,Oriol Vallès-Codina,         **
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
**           Authors: Diego A. Hartasánchez,Oriol Vallès-Codina,          **
**                   Marina Brasó-Vives, Juan Manuel Fuentes              **
**                   and Arcadi Navarro                                   **
**  Website/Contact: http://www.biologiaevolutiva.org/sedus/              **
**             Date: 01.10.14                                             **
**          Version: 1.10                                                 **
****************************************************************************/

#ifndef SEDUS_H
#define SEDUS_H

#include <QTextBrowser>
#include <QProgressBar>
#include <QMutex>
#include <QObject>
#include <QVector>
#include <stdbool.h>


typedef QVector<QVector<double>> qvdouble;

class sedus : public QObject
{

    Q_OBJECT
    struct prev_pres {
        int prev;
        int pres;
    };



public:
    struct param_exec{
        int runs;
        int sample_size;
        int snapshots;
    };
    struct param_main{
        int N;
        double theta;
        int blocklenght;
        int burnin;
        bool israndom;
        float fixation_linear;
        int total;
    };
    struct param_igc{
        float C;
        int lambda;
        float donor;
        int MEPS;
    };
    struct hotspot{
        double begin;
        double end;
        double rate;

    };

    struct param_crossover{
        float R;
        int hotspots_number;
        hotspot hotspots[5];
    };
    struct param_plot{
        bool piorS;
    };
    struct param_outs{
        bool proffile;
        bool pifile;
        bool Sfile;
        bool mutfile;
        bool SFSfile;
    };

    struct parameters{
        param_exec exec;
        param_main main;
        param_igc igc;
        param_crossover crossover;
        param_plot plot;
        param_outs outs;
        std::string id;
        std::string dir;
    };
    //extern parameters params;

    sedus(parameters *params, QObject *parent = 0);
    //sedus(QObject *parent = 0);
    //void dowork(int, char**, QTextBrowser*);
    void abort();

    void requestWork();
    void setLog(const QString &value);
    void setBar(int value);
    void setChart(const qvdouble &x,const qvdouble &y);
    bool _working;
private:

    bool _abort;
    QMutex mutex;
    void phaseIII();
    prev_pres phaseI();
    int phaseII(int,int,int);

signals:
    void workRequested();
    void valueChanged(const QString &value);
    void finished();
    void addLog(const QString &value);
    void addBar(int value);
    void addChart(const qvdouble &x,const qvdouble &y);
public slots:
    void dowork();


};

#endif // SEDUS_H


