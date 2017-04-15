#include "izhikevich_spnet.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    izhikevich_SPNET w;
    w.show();

    return a.exec();
}
