#ifndef HH_MODEL_QT_H
#define HH_MODEL_QT_H

#include <QMainWindow>
#include "qcustomplot.h"

namespace Ui {
class HH_Model_Qt;
}

class HH_Model_Qt : public QMainWindow
{
    Q_OBJECT

public:
    explicit HH_Model_Qt(QWidget *parent = 0);
    ~HH_Model_Qt();

    void Control(const int &cmd);
    void PrintPlot(QCustomPlot *HH_Plot);
    //void PrintPlot1(QCustomPlot *HH_Plot, QCustomPlot *HH_Plot1);

private slots:

    void on_horizontalSlider_sliderMoved();

    void on_horizontalSlider_2_sliderMoved();

    void on_horizontalSlider_3_sliderMoved();

    void on_horizontalSlider_4_sliderMoved();

    void on_horizontalSlider_5_sliderMoved();

    void on_pushButton_clicked();

private:
    Ui::HH_Model_Qt *ui;
    QString WinName;
};

#endif // HH_MODEL_QT_H
