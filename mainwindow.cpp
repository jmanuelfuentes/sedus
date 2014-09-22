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
    connect(ui->pushButton_dir, SIGNAL(released()), this, SLOT(dirButton()));
    connect(ui->fixlinearslider, SIGNAL(valueChanged(int)), this, SLOT(setFixLinearText(int)));
    connect(ui->fixlinearRadioButton, SIGNAL(toggled(bool)), this, SLOT(setFixLinearSliderVisible()));
    connect(ui->N, SIGNAL(textChanged(const QString)), this, SLOT(changeNdependecies(const QString)));

    //main window
    ui->id->setText("test-"+QDateTime::currentDateTime().toString("dd.MM.yyyy.hh.mm.ss"));
    ui->dir->setText(QDir::currentPath());
    ui->log->setVisible(false);
    ui->progressBar->setVisible(false);
    ui->customPlot->setVisible(false);
    ui->stop->setEnabled(false);

    //params execution
    ui->runs->setText("10");
    ui->samplesize->setText("50");
    ui->generations->setText("1000");
    //params main
    ui->N->setText("1000");
    ui->fixlinearRadioButton->setChecked(true);
    ui->theta->setText("0.001");
    ui->blocklength->setText("5000");
    ui->burnin->setText("30");
    ui->total->setText("130");

    //params igc
    ui->lambda->setText("100");
    ui->meps->setText("0");
    ui->C->setText("1");
    ui->donor->setText("0.5");

    //params crossover
    ui->R->setText("10");
    ui->hot1->setVisible(true);
    ui->hot2->setVisible(false);
    ui->hot3->setVisible(false);
    ui->hot4->setVisible(false);
    ui->hot5->setVisible(false);
    ui->hotspots->setValue(1);
    ui->rate1->setText("1");

    /*ui->R->setText("10");
    ui->C->setText("1");
    ui->superTime->setText("10");FET
    ui->timeToFixation->setText("20");FET
    ui->donoRatio->setText("0.5");FET
    ui->radioButton->setChecked(true);

*/

}

void MainWindow::hotspots(int value){
    switch (value)
    {
        case 1:{
            ui->hot1->setVisible(true);ui->hot2->setVisible(false);ui->hot3->setVisible(false);ui->hot4->setVisible(false);ui->hot5->setVisible(false);
            ui->rate1->setText("1");
            break;
        }
        case 2:{
            ui->hot1->setVisible(true);ui->hot2->setVisible(true);ui->hot3->setVisible(false);ui->hot4->setVisible(false);ui->hot5->setVisible(false);
            ui->rate1->setText("0.5");ui->rate2->setText("0.5");
            break;
        }
        case 3:{
            ui->hot1->setVisible(true);ui->hot2->setVisible(true);ui->hot3->setVisible(true);ui->hot4->setVisible(false);ui->hot5->setVisible(false);
            ui->rate1->setText("0.33");ui->rate2->setText("0.33");ui->rate3->setText("0.34");
            break;
        }
        case 4:{
            ui->hot1->setVisible(true);ui->hot2->setVisible(true);ui->hot3->setVisible(true);ui->hot4->setVisible(true);ui->hot5->setVisible(false);
            ui->rate1->setText("0.25");ui->rate2->setText("0.25");ui->rate3->setText("0.25");ui->rate4->setText("0.25");
            break;
        }
        case 5:{
            ui->hot1->setVisible(true);ui->hot2->setVisible(true);ui->hot3->setVisible(true);ui->hot4->setVisible(true);ui->hot5->setVisible(true);
            ui->rate1->setText("0.2");ui->rate2->setText("0.2");ui->rate3->setText("0.2");ui->rate4->setText("0.2");ui->rate5->setText("0.2");
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
    ui->customPlot->clearGraphs();
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
    ui->customPlot->clearGraphs();
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

bool MainWindow::checkRateHotSopts(){
    switch (ui->hotspots->value())
    {
        case 1:{
            sumratiohots=ui->rate1->text().toDouble();
            if(sumratiohots!=1) return false;
            break;
        }
        case 2:{
            sumratiohots=(ui->rate1->text().toDouble()+ui->rate2->text().toDouble());
            if(sumratiohots!=1) return false;
            break;
        }
        case 3:{
            sumratiohots=(ui->rate1->text().toDouble()+ui->rate2->text().toDouble()+ui->rate3->text().toDouble());
            if(sumratiohots!=1) return false;
            break;
        }
        case 4:{
            sumratiohots=(ui->rate1->text().toDouble()+ui->rate2->text().toDouble()+ui->rate3->text().toDouble()+ui->rate4->text().toDouble());
            std::cout << double(sumratiohots)<<" "<<double(1) << std::endl;
            if(sumratiohots!=1) return false;
            break;

        }
        case 5:{
            sumratiohots=(ui->rate1->text().toDouble()+ui->rate2->text().toDouble()+ui->rate3->text().toDouble()+ui->rate4->text().toDouble()+ui->rate4->text().toDouble());
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
    if(!checkRateHotSopts()){
        QMessageBox::information(this, tr("ERROR HOTSPOTS"), "The sum of  the all hotspots ratio must be one. (Current sum: "+QString::number(sumratiohots)+")");
        return;
    }

    if(!ui->progressBar->isVisible()){
        ui->progressBar->setMaximum(100);
        ui->progressBar->setValue(0);
        ui->progressBar->setVisible(true);
    }
    ui->customPlot->setVisible(true);
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
    params.exec.sample_size = ui->runs->text().toInt();
    params.exec.generations = ui->generations->text().toInt();
    //params main
    params.main.N = ui->N->text().toInt();
    params.main.theta = ui->theta->text().toDouble();
    params.main.blocklenght = ui->blocklength->text().toInt();
    params.main.burnin = ui->burnin->text().toInt();
    params.main.total = ui->total->text().toInt();
    if(ui->fixlinearRadioButton->isChecked()){
        params.main.israndom = false;
        params.main.fixation_linear = ui->fixlinearlabel->text().toFloat();
    }else{
        params.main.israndom = true;
    }
    //params igc
    params.igc.C = ui->C->text().toFloat();
    params.igc.donor = ui->donor->text().toFloat();
    params.igc.lambda = ui->lambda->text().toInt();
    params.igc.MEPS = ui->meps->text().toInt();
    //params crossover
    params.crossover.R = ui->R->text().toFloat();



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
    ui->customPlot->graph(0)->setPen(QPen(Qt::darkGreen));
    //ui->customPlot->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20)));
    ui->customPlot->addGraph();
    ui->customPlot->graph(1)->setName("Single-copy");
    ui->customPlot->graph(1)->setPen(QPen(Qt::red));
    ui->customPlot->addGraph();
    ui->customPlot->graph(2)->setName("Duplicated");
    ui->customPlot->graph(2)->setPen(QPen(QColor(1,0,1)));
    ui->customPlot->addGraph();
    ui->customPlot->graph(3)->setName("O+D");
    ui->customPlot->graph(3)->setPen(QPen(Qt::blue));

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

    ui->customPlot->xAxis->setRange(0, (x[0].length()*ui->generations->text().toInt())/1000);
    ui->customPlot->yAxis->setRange(0, 15);
    ui->customPlot->graph(0)->setData(x[0], y[0]);
    ui->customPlot->graph(1)->setData(x[1], y[1]);
    ui->customPlot->graph(2)->setData(x[2], y[2]);
    ui->customPlot->graph(3)->setData(x[3], y[3]);
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
}
