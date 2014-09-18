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

private:
    QThread *thread;
    sedus *worker;

    //QPushButton *m_button;
};

#endif // MAINWINDOW_H
