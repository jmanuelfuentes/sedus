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
#include "ui_mainwindow.h"
#include <QDesktopWidget>
#include <QApplication>
#include <qthread.h>
#include <sedus.h>
#include <iostream>
#include <sstream>
#include <string>
#include <QFileDialog>
#include <QDateTime>
#include <stdbool.h>

const char kPathSeparator =
#ifdef _WIN32
                            '\\';
#else
                            '/';
#endif



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    uidelete = false;
    ui->setupUi(this);
    // Create the button, make "this" the parent
        //m_button = new QPushButton("My Button", this);
        // set size and location of the button
        //m_button->setGeometry(QRect(QPoint(230, 180),QSize(99, 27)));

        // Connect button signal to appropriate slot
    centerWidget();
    connect(ui->m_button, SIGNAL(released()), this, SLOT(handleButton()));
    connect(ui->stop, SIGNAL(released()), this, SLOT(stop()));
    connect(ui->hotspots, SIGNAL(valueChanged(int)), this, SLOT(hotspots(int)));
    connect(ui->R, SIGNAL(textChanged(const QString)), this, SLOT(changeR(const QString)));
    connect(ui->pushButton_dir, SIGNAL(released()), this, SLOT(dirButton()));
    connect(ui->fixlinearslider, SIGNAL(valueChanged(int)), this, SLOT(setFixLinearText(int)));
    connect(ui->burnin_slider, SIGNAL(valueChanged(int)), this, SLOT(setBurnin(int)));
    connect(ui->total_slider, SIGNAL(valueChanged(int)), this, SLOT(setTotal(int)));
    connect(ui->fixlinearRadioButton, SIGNAL(toggled(bool)), this, SLOT(setFixLinearSliderVisible()));
    connect(ui->N, SIGNAL(textChanged(const QString)), this, SLOT(changeNdependecies(const QString)));
    connect(ui->snapshots, SIGNAL(textChanged(const QString)), this, SLOT(changeGenerationsDependecies(const QString)));
    connect(this, SIGNAL(setMinimum(double,int)), this, SLOT(changeMinimum(double,int)));
    connect(this, SIGNAL(setMaximum(double,int)), this, SLOT(changeMaximum(double,int)));
    connect(ui->b1, SIGNAL(valueChanged(double)), this, SLOT(setMinimumE1(double)));
    connect(ui->e1, SIGNAL(valueChanged(double)), this, SLOT(setMaximumE1(double)));
    connect(ui->b2, SIGNAL(valueChanged(double)), this, SLOT(setMinimumE2(double)));
    connect(ui->e2, SIGNAL(valueChanged(double)), this, SLOT(setMaximumE2(double)));
    connect(ui->b3, SIGNAL(valueChanged(double)), this, SLOT(setMinimumE3(double)));
    connect(ui->e3, SIGNAL(valueChanged(double)), this, SLOT(setMaximumE3(double)));
    connect(ui->b4, SIGNAL(valueChanged(double)), this, SLOT(setMinimumE4(double)));
    connect(ui->e4, SIGNAL(valueChanged(double)), this, SLOT(setMaximumE4(double)));
    connect(ui->b5, SIGNAL(valueChanged(double)), this, SLOT(setMinimumE5(double)));
    connect(ui->e5, SIGNAL(valueChanged(double)), this, SLOT(setMaximumE5(double)));
    connect(ui->rate1, SIGNAL(valueChanged(double)), this, SLOT(changeAlpha1(double)));
    connect(ui->rate2, SIGNAL(valueChanged(double)), this, SLOT(changeAlpha2(double)));
    connect(ui->rate3, SIGNAL(valueChanged(double)), this, SLOT(changeAlpha3(double)));
    connect(ui->rate4, SIGNAL(valueChanged(double)), this, SLOT(changeAlpha4(double)));
    connect(ui->rate5, SIGNAL(valueChanged(double)), this, SLOT(changeAlpha5(double)));

    //main window
    ui->id->setText("test-"+QDateTime::currentDateTime().toString("dd.MM.yyyy"));
    ui->dir->setText(QDir::currentPath());
    ui->log->setVisible(false);
    ui->progressBar->setVisible(false);
    ui->customPlot->setVisible(false);
    ui->stop->setEnabled(false);

    //params execution
    ui->runs->setText("10");
    ui->samplesize->setText("50");
    ui->snapshots->setText("1000");

    //params main
    ui->N->setText("1000");
    ui->radioButton->setChecked(true);
    ui->fixgroupslider->setVisible(false);
    ui->theta->setText("0.001");
    ui->blocklength->setText("5000");

    //params igc
    ui->lambda->setText("100");
    ui->meps->setText("0");
    ui->C->setText("1");
    //ui->donor->setText("0.5");
    ui->donor_spin->setValue(0.5);

    //params crossover
    ui->R->setText("1");
    ui->hot1->setVisible(true);
    ui->hot2->setVisible(false);
    ui->hot3->setVisible(false);
    ui->hot4->setVisible(false);
    ui->hot5->setVisible(false);
    ui->lh1->setVisible(true);
    ui->lh2->setVisible(false);
    ui->lh3->setVisible(false);
    ui->lh4->setVisible(false);
    ui->lh5->setVisible(false);
    ui->hotspots->setValue(1);
    ui->rate1->setValue(1.00);

    //params plot
    ui->plotpi->setChecked(true);

    //params outs
    ui->proffile->setChecked(true);
    ui->pifile->setChecked(true);
    ui->Sfile->setChecked(true);
    ui->mutfile->setChecked(true);
    ui->SFSfile->setChecked(true);


 //   QPalette palette = ui->colh1->palette();
 //   QColor color = palette.color(QPalette::Base);
 //   color.setAlpha(255); // 1 to 255
 //   palette.setColor(QPalette::Base, color);
 //   ui->colh1->setPalette(palette);


    /*ui->R->setText("10");
    ui->C->setText("1");
    ui->superTime->setText("10");FET
    ui->timeToFixation->setText("20");FET
    ui->donoRatio->setText("0.5");FET
    ui->radioButton->setChecked(true);

*/

}

void MainWindow::setMinimumE1(double value){
    emit setMinimum(value,1);
}
void MainWindow::setMaximumE1(double value){
    emit setMaximum(value,1);
}
void MainWindow::setMinimumE2(double value){
    emit setMinimum(value,2);
}
void MainWindow::setMaximumE2(double value){
    emit setMaximum(value,2);
}
void MainWindow::setMinimumE3(double value){
    emit setMinimum(value,3);
}
void MainWindow::setMaximumE3(double value){
    emit setMaximum(value,3);
}
void MainWindow::setMinimumE4(double value){
    emit setMinimum(value,4);
}
void MainWindow::setMaximumE4(double value){
    emit setMaximum(value,4);
}
void MainWindow::setMinimumE5(double value){
    emit setMinimum(value,5);
}
void MainWindow::setMaximumE5(double value){
    emit setMaximum(value,5);
}

void MainWindow::changeAlpha1(double value){
       QPalette palette = ui->ch1->palette();
       QColor color = palette.color(QPalette::Base);
       color.setAlpha((int)(value*255)); // 1 to 255
       palette.setColor(QPalette::Base, color);
       ui->ch1->setPalette(palette);
       ui->lh1->setPalette(palette);
}
void MainWindow::changeAlpha2(double value){
       QPalette palette = ui->ch2->palette();
       QColor color = palette.color(QPalette::Base);
       color.setAlpha((int)(value*255)); // 1 to 255
       palette.setColor(QPalette::Base, color);
       ui->ch2->setPalette(palette);
       ui->lh2->setPalette(palette);
}
void MainWindow::changeAlpha3(double value){
       QPalette palette = ui->ch3->palette();
       QColor color = palette.color(QPalette::Base);
       color.setAlpha((int)(value*255)); // 1 to 255
       palette.setColor(QPalette::Base, color);
       ui->ch3->setPalette(palette);
       ui->lh3->setPalette(palette);
}
void MainWindow::changeAlpha4(double value){
       QPalette palette = ui->ch4->palette();
       QColor color = palette.color(QPalette::Base);
       color.setAlpha((int)(value*255)); // 1 to 255
       palette.setColor(QPalette::Base, color);
       ui->ch4->setPalette(palette);
       ui->lh4->setPalette(palette);
}
void MainWindow::changeAlpha5(double value){
       QPalette palette = ui->ch5->palette();
       QColor color = palette.color(QPalette::Base);
       color.setAlpha((int)(value*255)); // 1 to 255
       palette.setColor(QPalette::Base, color);
       ui->ch5->setPalette(palette);
       ui->lh5->setPalette(palette);
}


void MainWindow::changeMinimum(double value, int item){
    switch (item){
        case 1:{
            if(value<3)ui->e1->setMinimum(value+0.0001);
            else ui->b1->setValue(ui->e1->value()-0.0001);
            ui->lh1->setFixedWidth((ui->e1->value() - ui->b1->value())*100);
            ui->lh1->setGeometry(ui->b1->value()*100+ui->blocks->geometry().x(),ui->lh1->geometry().y(),ui->lh1->geometry().width(),ui->lh1->geometry().height());
            break;
        }
    case 2:{
        if(value<3)ui->e2->setMinimum(value+0.0001);
        else ui->b2->setValue(ui->e2->value()-0.0001);
        ui->lh2->setFixedWidth((ui->e2->value() - ui->b2->value())*100);
        ui->lh2->setGeometry(ui->b2->value()*100+ui->blocks->geometry().x(),ui->lh2->geometry().y(),ui->lh2->geometry().width(),ui->lh2->geometry().height());
        break;
    }
    case 3:{
        if(value<3)ui->e3->setMinimum(value+0.0001);
        else ui->b3->setValue(ui->e3->value()-0.0001);
        ui->lh3->setFixedWidth((ui->e3->value() - ui->b3->value())*100);
        ui->lh3->setGeometry(ui->b3->value()*100+ui->blocks->geometry().x(),ui->lh3->geometry().y(),ui->lh3->geometry().width(),ui->lh3->geometry().height());
        break;
    }
    case 4:{
       if(value<3) ui->e4->setMinimum(value+0.0001);
       else ui->b4->setValue(ui->e4->value()-0.0001);
       ui->lh4->setFixedWidth((ui->e4->value() - ui->b4->value())*100);
       ui->lh4->setGeometry(ui->b4->value()*100+ui->blocks->geometry().x(),ui->lh4->geometry().y(),ui->lh4->geometry().width(),ui->lh4->geometry().height());
        break;
    }
    case 5:{
        if(value<3)ui->e5->setMinimum(value+0.0001);
        else ui->b5->setValue(ui->e5->value()-0.0001);
        ui->lh5->setFixedWidth((ui->e5->value() - ui->b5->value())*100);
        ui->lh5->setGeometry(ui->b5->value()*100+ui->blocks->geometry().x(),ui->lh5->geometry().y(),ui->lh5->geometry().width(),ui->lh5->geometry().height());
        break;
    }

    }
}
void MainWindow::changeMaximum(double value, int item){
    switch (item){
    case 1:{
            if(value>0)ui->b1->setMaximum(value-0.0001);
            else ui->e1->setValue(ui->b1->value()+0.0001);
            ui->lh1->setVisible(false);
            ui->lh1->setVisible(true);
            ui->lh1->setFixedWidth((ui->e1->value() - ui->b1->value())*100);
            ui->lh1->setGeometry(ui->b1->value()*100+ui->blocks->geometry().x(),ui->lh1->geometry().y(),ui->lh1->geometry().width(),ui->lh1->geometry().height());
            break;
        }
    case 2:{
        if(value>0)ui->b2->setMaximum(value-0.0001);
        else ui->e2->setValue(ui->b2->value()+0.0001);
        ui->lh2->setVisible(false);
        ui->lh2->setVisible(true);
        ui->lh2->setFixedWidth((ui->e2->value() - ui->b2->value())*100);
        break;
    }
    case 3:{
        if(value>0)ui->b3->setMaximum(value-0.0001);
        else ui->e3->setValue(ui->b3->value()+0.0001);
        ui->lh3->setVisible(false);
        ui->lh3->setVisible(true);
        ui->lh3->setFixedWidth((ui->e3->value() - ui->b3->value())*100);
        break;
    }
    case 4:{
        if(value>0)ui->b4->setMaximum(value-0.0001);
        else ui->e4->setValue(ui->b4->value()+0.0001);
        ui->lh4->setVisible(false);
        ui->lh4->setVisible(true);
        ui->lh4->setFixedWidth((ui->e4->value() - ui->b4->value())*100);
        break;
    }
    case 5:{
        if(value>0)ui->b5->setMaximum(value-0.0001);
        else ui->e5->setValue(ui->b5->value()+0.0001);
        ui->lh5->setVisible(false);
        ui->lh5->setVisible(true);
        ui->lh5->setFixedWidth((ui->e5->value() - ui->b5->value())*100);
        break;
    }
    }
}

void MainWindow::hotspots(int value){
    switch (value)
    {
        case 1:{
            ui->hot1->setVisible(true);ui->hot2->setVisible(false);ui->hot3->setVisible(false);ui->hot4->setVisible(false);ui->hot5->setVisible(false);
            ui->rate1->setValue(1);
            ui->lh1->setVisible(true);ui->lh2->setVisible(false);ui->lh3->setVisible(false);ui->lh4->setVisible(false);ui->lh5->setVisible(false);

            break;
        }
        case 2:{
            ui->hot1->setVisible(true);ui->hot2->setVisible(true);ui->hot3->setVisible(false);ui->hot4->setVisible(false);ui->hot5->setVisible(false);
            ui->rate1->setValue(0.5);ui->rate2->setValue(0.5);
            ui->lh1->setVisible(true);ui->lh2->setVisible(true);ui->lh3->setVisible(false);ui->lh4->setVisible(false);ui->lh5->setVisible(false);

            break;
        }
        case 3:{
            ui->hot1->setVisible(true);ui->hot2->setVisible(true);ui->hot3->setVisible(true);ui->hot4->setVisible(false);ui->hot5->setVisible(false);
            ui->rate1->setValue(0.33);ui->rate2->setValue(0.33);ui->rate3->setValue(0.34);
            ui->lh1->setVisible(true);ui->lh2->setVisible(true);ui->lh3->setVisible(true);ui->lh4->setVisible(false);ui->lh5->setVisible(false);

            break;
        }
        case 4:{
            ui->hot1->setVisible(true);ui->hot2->setVisible(true);ui->hot3->setVisible(true);ui->hot4->setVisible(true);ui->hot5->setVisible(false);
            ui->rate1->setValue(0.25);ui->rate2->setValue(0.25);ui->rate3->setValue(0.25);ui->rate4->setValue(0.25);
            ui->lh1->setVisible(true);ui->lh2->setVisible(true);ui->lh3->setVisible(true);ui->lh4->setVisible(true);ui->lh5->setVisible(false);

            break;
        }
        case 5:{
            ui->hot1->setVisible(true);ui->hot2->setVisible(true);ui->hot3->setVisible(true);ui->hot4->setVisible(true);ui->hot5->setVisible(true);
            ui->rate1->setValue(0.2);ui->rate2->setValue(0.2);ui->rate3->setValue(0.2);ui->rate4->setValue(0.2);ui->rate5->setValue(0.2);
            ui->lh1->setVisible(true);ui->lh2->setVisible(true);ui->lh3->setVisible(true);ui->lh4->setVisible(true);ui->lh5->setVisible(true);

            break;
        }
    }
}

void MainWindow::stop(){
    //worker->abort();
    //delete worker;
    worker->abort();
    //thread->quit();


    ui->progressBar->setValue(0);
    ui->m_button->setEnabled(true);
    ui->stop->setEnabled(false);
    ui->log->clear();
    //ui->customPlot->clearGraphs();
    //ui->customPlot->setEnabled(false);
}

void MainWindow::quitWork(){
    thread->quit();
    thread->wait();
    delete thread;
    delete worker;

    ui->progressBar->setValue(0);
    ui->m_button->setEnabled(true);
    ui->stop->setEnabled(false);
    ui->log->clear();
    //ui->customPlot->clearGraphs();
    //ui->customPlot->setEnabled(false);
    if(uidelete){delete ui;}
}

MainWindow::~MainWindow()
{
    if(thread != NULL){
        worker->abort();
        uidelete = true;
    }else{
        delete ui;
    }

}
void MainWindow::dirButton()
{
    QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),"/home",QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
   ui->dir->setText(dir);
}

bool MainWindow::checkN(){

    if(!(ui->N->text().toInt()%2==0))return false;
    return true;
}
bool MainWindow::checkGenerations(){

    if(!(ui->snapshots->text().toInt()<=ui->N->text().toInt()))return false;
    if(!(ui->snapshots->text().toInt()%2==0))return false;
    if(!(ui->N->text().toInt()%ui->snapshots->text().toInt())==0)return false;
    return true;
}

bool MainWindow::checkRateHotSopts(){
    switch (ui->hotspots->value())
    {
        case 1:{
            sumratiohots=ui->rate1->value();
            if(sumratiohots!=1) return false;
            break;
        }
        case 2:{
            sumratiohots=(ui->rate1->value()+ui->rate2->value());
            if(sumratiohots!=1) return false;
            break;
        }
        case 3:{
            sumratiohots=(ui->rate1->value()+ui->rate2->value()+ui->rate3->value());
            if(sumratiohots!=1) return false;
            break;
        }
        case 4:{
            sumratiohots=(ui->rate1->value()+ui->rate2->value()+ui->rate3->value()+ui->rate4->value());
            std::cout << double(sumratiohots)<<" "<<double(1) << std::endl;
            if(sumratiohots!=1) return false;
            break;

        }
        case 5:{
            sumratiohots=(ui->rate1->value()+ui->rate2->value()+ui->rate3->value()+ui->rate4->value()+ui->rate4->value());
            if(sumratiohots!=1) return false;
            break;
        }

     }
    return true;
}

void MainWindow::handleButton()
{
    typedef QVector<QVector<double>> qvdouble;

    qRegisterMetaType<qvdouble>("qvdouble");

    //CHECKS
    if(!checkN()){
        QMessageBox::information(this, tr("ERROR N"), "\"Population size (N)\" must be even.");
        return;
    }
    if(!checkGenerations()){
        QMessageBox::information(this, tr("ERROR GENERATIONS"), "\"Generations between snapshots\" must be even, less than or equal to \"Population size (N)\" and a factor of N.");
        return;
    }
    if(!checkRateHotSopts()){
        QMessageBox::information(this, tr("ERROR HOTSPOTS"), "The sum of the rations for all crossover regions must be 1. (Current sum: "+QString::number(sumratiohots)+")");
        return;
    }


    if(!ui->progressBar->isVisible()){
        ui->progressBar->setMaximum(100);
        ui->progressBar->setValue(0);
        ui->progressBar->setVisible(true);
    }
    ui->customPlot->setVisible(true);
    ui->customPlot->clearGraphs();
    ui->log->setVisible(true);
    ui->m_button->setEnabled(false);
    ui->stop->setEnabled(true);


    //argv[1]=ui->letter->text().toLocal8Bit().data();
    //QString s = QString::fromStdString(argv[1]);
    //params = new parameters;

    sedus::parameters params;
    //main window
    params.id = ui->id->text().toLocal8Bit().data();
    params.dir = ui->dir->text().toLocal8Bit().data();
    //params exec
    params.exec.runs = ui->runs->text().toInt();
    params.exec.sample_size = ui->samplesize->text().toInt();
    params.exec.snapshots = ui->snapshots->text().toInt();
    //params main
    params.main.N = ui->N->text().toInt();
    params.main.theta = ui->theta->text().toDouble();
    params.main.blocklenght = ui->blocklength->text().toInt();
    params.main.burnin = ui->burnin_slider->value();
    params.main.total = ui->total_slider->value();
    if(ui->fixlinearRadioButton->isChecked()){
        params.main.israndom = false;
        params.main.fixation_linear = ui->fixlinearlabel->text().toFloat();
    }else{
        params.main.israndom = true;
    }
    if(ui->proffile->isChecked()){params.outs.proffile = true;}
    if(ui->pifile->isChecked()){params.outs.pifile = true;}
    if(ui->Sfile->isChecked()){params.outs.Sfile = true;}
    if(ui->mutfile->isChecked()){params.outs.mutfile = true;}
    if(ui->SFSfile->isChecked()){params.outs.SFSfile = true;}

    //params igc
    params.igc.C = ui->C->text().toFloat();
    //params.igc.donor = ui->donor->text().toFloat();
    params.igc.donor = ui->donor_spin->value();
    params.igc.lambda = ui->lambda->text().toInt();
    params.igc.MEPS = ui->meps->text().toInt();
    //params crossover
    params.crossover.R = ui->R->text().toFloat();
    if(params.crossover.R!=float(0)){
        params.crossover.hotspots_number = ui->hotspots->value();
        for(int i=1; i<=params.crossover.hotspots_number;i++){
            switch (i)
            {
                case 1:
                params.crossover.hotspots[0].begin = ui->b1->value();
                params.crossover.hotspots[0].end = ui->e1->value();
                params.crossover.hotspots[0].rate = ui->rate1->value();
                break;
            case 2:
            params.crossover.hotspots[1].begin = ui->b2->value();
            params.crossover.hotspots[1].end = ui->e2->value();
            params.crossover.hotspots[1].rate = ui->rate2->value();
            break;
            case 3:
            params.crossover.hotspots[2].begin = ui->b3->value();
            params.crossover.hotspots[2].end = ui->e3->value();
            params.crossover.hotspots[2].rate = ui->rate3->value();
            break;
            case 4:
            params.crossover.hotspots[3].begin = ui->b4->value();
            params.crossover.hotspots[3].end = ui->e4->value();
            params.crossover.hotspots[3].rate = ui->rate4->value();
            break;
            case 5:
            params.crossover.hotspots[4].begin = ui->b5->value();
            params.crossover.hotspots[4].end = ui->e5->value();
            params.crossover.hotspots[4].rate = ui->rate5->value();
            break;
            }

        }
    }
    //params plot
    if(ui->plotpi->isChecked()){
        params.plot.piorS = true;
    }else{
        params.plot.piorS = false;
    }


    #ifdef _WIN32
        // replace(argv[7], "\\", "\\\\");
        params.dir=params.dir+"\\";
    #else
        params.dir=params.dir+"/";
    #endif

    /*if(ui->radioButton->isChecked()){
        argv[8]="SC";
    }else if(ui->radioButton_2->isChecked()){
        argv[8]="WR";
    }else if(ui->radioButton_3->isChecked()){
        argv[8]="HS";
    }*/

    //chart
    ui->customPlot->legend->setVisible(true);
    ui->customPlot->legend->setFont(QFont("Helvetica",9));
    ui->customPlot->axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignLeft|Qt::AlignTop);
    ui->customPlot->addGraph();

    ui->customPlot->graph(0)->setName("Original");
    ui->customPlot->graph(0)->setPen(QPen(QColor(26,30,178)));
    ui->customPlot->addGraph();
    ui->customPlot->graph(1)->setName("Single-copy");
    ui->customPlot->graph(1)->setPen(QPen(QColor(255,124,0)));
    ui->customPlot->addGraph();
    ui->customPlot->graph(2)->setName("Duplicated");
    ui->customPlot->graph(2)->setPen(QPen(QColor(62,148,209)));

    /*ui->customPlot->addGraph();
    ui->customPlot->graph(3)->setName("O+D");
    ui->customPlot->graph(3)->setPen(QPen(Qt::blue));
    */
    ui->customPlot->xAxis2->setVisible(true);
    ui->customPlot->xAxis2->setTickLabels(false);
    ui->customPlot->yAxis2->setVisible(true);
    ui->customPlot->yAxis2->setTickLabels(false);
      // make left and bottom axes always transfer their ranges to right and top axes:
    connect(ui->customPlot->xAxis, SIGNAL(rangeChanged(QCPRange)), ui->customPlot->xAxis2, SLOT(setRange(QCPRange)));
    connect(ui->customPlot->yAxis, SIGNAL(rangeChanged(QCPRange)), ui->customPlot->yAxis2, SLOT(setRange(QCPRange)));
      // pass data points to graphs:
    //ui->customPlot->graph(0)->setData(x, y0);
      // let the ranges scale themselves so graph 0 fits perfectly in the visible area:
    ui->customPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    //ui->customPlot->replot();

    thread = new QThread();
    worker = new sedus(&params);
    //worker = new sedus();
    worker->moveToThread(thread);
    connect(worker, SIGNAL(valueChanged(QString)), ui->label, SLOT(setText(QString)));
    connect(worker, SIGNAL(workRequested()), thread, SLOT(start()));
    connect(thread, SIGNAL(started()), worker, SLOT(dowork()));
    connect(worker, SIGNAL(finished()), this, SLOT(quitWork()));
    //connect(worker, SIGNAL(finished()), thread, SLOT(quit()), Qt::DirectConnection);
    connect(worker, SIGNAL(addLog(const QString)),this, SLOT(setLog(const QString)));
    connect(worker, SIGNAL(addBar(int)),this, SLOT(setBar(int)));
    connect(worker, SIGNAL(addChart(qvdouble,qvdouble)),this, SLOT(setChart(qvdouble,qvdouble)));

    //sedus(8,argv,ui->log);
    worker->abort();
    thread->wait(); // If the thread is not running, this will immediately return.
    worker->requestWork();

    /*for(int i=0; i<100; i++){
        ui->progressBar->setValue(i);
        QThread::msleep(1000);
    }*/
}

void MainWindow::setChart(const qvdouble &x,const qvdouble &y){

    ui->customPlot->xAxis->setRange(0, (x[0].length()*ui->snapshots->text().toInt())/1000);
    ui->customPlot->yAxis->setRange(0, 20);
    ui->customPlot->xAxis->setLabel("Thousands of generations");
    ui->customPlot->yAxis->setLabel("Average pairwise differences");
    if(ui->plotpi->isChecked()){
        ui->customPlot->graph(0)->setData(x[0], y[0]);
        ui->customPlot->graph(1)->setData(x[1], y[1]);
        ui->customPlot->graph(2)->setData(x[2], y[2]);
    }
    if(ui->plotS->isChecked()){
        ui->customPlot->graph(0)->setData(x[0], y[3]);
        ui->customPlot->graph(1)->setData(x[1], y[4]);
        ui->customPlot->graph(2)->setData(x[2], y[5]);
    }
    //ui->customPlot->graph(3)->setData(x[3], y[3]);
    ui->customPlot->replot();
}
void MainWindow::setLog(const QString &value){
    ui->log->append(value);
}
void MainWindow::setBar(int value){
    ui->progressBar->setValue(value);
}
void MainWindow::centerWidget()
{
    QDesktopWidget* desk = QApplication::desktop();
    QRect screenGeometry = desk->screenGeometry();
    int x = (screenGeometry.width()-this->width()) / 2;
    int y = (screenGeometry.height()-this->height()) / 2;
    this->move(x, y);
    this->show();
}

bool replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

void MainWindow::setFixLinearText(int value){
    ui->fixlinearlabel->setText(QString::number(value));
}

void MainWindow::setBurnin(int value){
    ui->burnin_label->setText(QString::number(value));
    ui->total_slider->setMinimum(ui->burnin_slider->value()+20*ui->N->text().toInt());
    ui->total_slider->setMaximum(ui->N->text().toInt()*1000);
}

void MainWindow::setTotal(int value){
    ui->total_label_2->setText(QString::number(value));
}

void MainWindow::setFixLinearSliderVisible(){
    if (ui->fixlinearRadioButton->isChecked()){
        ui->fixlinearslider->setMaximum(20*(ui->N->text().toInt()));
        ui->fixlinearslider->setMinimum(1);
        ui->fixlinearslider->setValue(1*(ui->N->text().toInt()));
        ui->fixgroupslider->setVisible(true);
    }else{
        ui->fixgroupslider->setVisible(false);
    }
}
void MainWindow::changeNdependecies(const QString &value){
    ui->fixlinearslider->setMaximum(20*(ui->N->text().toInt()));
    ui->fixlinearslider->setMinimum(1);
    ui->fixlinearslider->setValue(1*(value.toInt()));

    //ui->burnin_slider->setMaximum(50*(ui->N->text().toInt()));
    ui->burnin_slider->setMaximum((50*(ui->N->text().toInt()))-(50*(ui->N->text().toInt())%ui->snapshots->text().toInt()));
    ui->burnin_slider->setMinimum(ui->snapshots->text().toInt());
    ui->burnin_slider->setValue(30*(value.toInt()));
    ui->burnin_slider->setSingleStep(ui->snapshots->text().toInt());

    //ui->total_slider->setMaximum(1000*(ui->N->text().toInt()));
    ui->total_slider->setMaximum((1000*(ui->N->text().toInt()))-(50*(ui->N->text().toInt())%ui->snapshots->text().toInt()));
    ui->total_slider->setMinimum(ui->burnin_slider->value()+20*ui->N->text().toInt());
    ui->total_slider->setValue(130*(value.toInt()));
    ui->total_slider->setSingleStep(ui->snapshots->text().toInt());
}
void MainWindow::changeGenerationsDependecies(const QString &value){
    if(value.toInt()>0){
        ui->burnin_slider->setMaximum((50*(ui->N->text().toInt()))-(50*(ui->N->text().toInt())%ui->snapshots->text().toInt()));
        ui->burnin_slider->setMinimum(value.toInt());
        ui->burnin_slider->setSingleStep(value.toInt());
        ui->burnin_slider->setValue((30*(ui->N->text().toInt()))-(30*(ui->N->text().toInt())%ui->snapshots->text().toInt()));
        ui->total_slider->setMaximum((1000*(ui->N->text().toInt()))-(1000*(ui->N->text().toInt())%ui->snapshots->text().toInt()));
        ui->total_slider->setSingleStep(value.toInt());
        ui->total_slider->setValue((130*(ui->N->text().toInt()))-(130*(ui->N->text().toInt())%ui->snapshots->text().toInt()));
    }else{
        //ui->snapshots->setText("1000");
    }
}

void MainWindow::changeR(const QString &value){
    /*if(value.compare("")==0){
        ui->R->setText("0");
    }*/
    if(value.toInt()!=0){
        ui->frame_2->setVisible(true);
        ui->groupBox_2->setVisible(true);
    }else{
        ui->frame_2->setVisible(false);
        ui->groupBox_2->setVisible(false);
    }

}
