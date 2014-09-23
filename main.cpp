#include "mainwindow.h"
#include <QApplication>
#include <QSplashScreen>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    QPixmap pixmap("./img/splash.png");
    QSplashScreen splash(pixmap);
    splash.show();
    splash.showMessage(QObject::tr("Initiating SEDUS..."),
                        Qt::AlignLeft | Qt::AlignTop, Qt::black);

    a.processEvents();
    QThread::sleep(5);
    MainWindow w;
    w.show();
    splash.finish(&w);
    splash.raise();
    return a.exec();
}
