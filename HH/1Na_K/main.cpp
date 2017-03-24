#include "hh_model_qt.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    HH_Model_Qt w;
    w.show();

    return a.exec();
}
