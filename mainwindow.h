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
    void dirButton();
    void centerWidget();
    void setLog(const QString &value);
    void setBar(int value);
    void setChart(const qvdouble &x,const qvdouble &y);
    void setFixLinearText(int value);
    void setFixLinearSliderVisible();
    void changeNdependecies(const QString &value);
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

signals:
    void setMinimum(double value, int item);
    void setMaximum(double value, int item);
private:
    QThread *thread;
    sedus *worker;
    bool checkRateHotSopts();
    bool uidelete;
    double sumratiohots;
    //QPushButton *m_button;
};

#endif // MAINWINDOW_H
