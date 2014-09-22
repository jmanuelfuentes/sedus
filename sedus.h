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
        int generations;
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
        double start;
        double end;
        double rate;

    };

    struct param_crossover{
        float R;
        bool isSC;
        bool isWR;
        bool isHS;
        int hotspots_number;
        hotspot hotspots[5];
    };

    struct parameters{
        param_exec exec;
        param_main main;
        param_igc igc;
        param_crossover crossover;
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


