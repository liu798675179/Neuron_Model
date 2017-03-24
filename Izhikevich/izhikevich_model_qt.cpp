#include "izhikevich_model_qt.h"
#include "ui_izhikevich_model_qt.h"
#include "Neuron.h"
#include <stdexcept>

using std::runtime_error;

HH_Model_Qt::HH_Model_Qt(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::HH_Model_Qt)
{
    ui->setupUi(this);
    Control(0);
}

HH_Model_Qt::~HH_Model_Qt()
{
    delete ui;
}

void HH_Model_Qt::Control(const int &cmd){
    if(cmd == 0){
        ui->HH_Plot->clearGraphs();
        PrintPlot(ui->HH_Plot);
        setWindowTitle(WinName);
        ui->HH_Plot->replot();
    }
    else if(cmd == 1){
        ui->HH_Plot->clearGraphs();
        setWindowTitle(WinName);
        ui->HH_Plot->replot();
    }
    else{
        try{
            throw runtime_error("Commend error.");
        }
        catch(runtime_error const  &err){
            err.what();
        }
    }
}

auto num = T / 8 / dt;

void HH_Model_Qt::PrintPlot(QCustomPlot *HH_Plot){
    WinName = "Izhikevich model.";

    auto t = 0.0;
    QVector<double> x(num), y0(num);

    Neuron temp_Neuron;
    auto get_Neuron = temp_Neuron.Simulation();

    // add two new graphs and set their look:
    HH_Plot->addGraph(HH_Plot->xAxis, HH_Plot->yAxis);
    HH_Plot->graph(0)->setPen(QPen(Qt::black));
    HH_Plot->graph(0)->setName("The membrane potential.");
    HH_Plot->legend->setVisible(true);
    HH_Plot->legend->setBrush(QColor(0,0,0,0));

    // generate some points of data (y0 for first, y1 for second graph):
    for(auto i = 0; i < num; ++i){
        x[i] = t;
        y0[i] = get_Neuron[i];

        t += dt;
    }

    // configure right and top axis to show ticks but no labels:
    HH_Plot->xAxis->setLabel("Time: mS");
    HH_Plot->xAxis2->setLabel("Izhikevich");
    HH_Plot->yAxis->setLabel("Potential: mV");


    // make left and bottom axes always transfer their ranges to right and top axes:
    connect(HH_Plot->xAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot->xAxis2, SLOT(setRange(QCPRange)));
    connect(HH_Plot->yAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot->yAxis2, SLOT(setRange(QCPRange)));

    // pass data points to graphs:
    HH_Plot->graph(0)->addData(x, y0);


    // let the ranges scale themselves so graph 0 fits perfectly in the visible area:
    HH_Plot->graph(0)->rescaleAxes();

    // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
    HH_Plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
}

void HH_Model_Qt::on_horizontalSlider_sliderMoved() {
    ui->horizontalSlider->setRange(5, T);
    ui->label_time->setText(QString::number(ui->horizontalSlider->value()) + " ms");
    num = ui->horizontalSlider->value() / dt;
    Control(0);
}

void HH_Model_Qt::on_horizontalSlider_2_sliderMoved() {
    ui->horizontalSlider_2->setRange(-10, 110);

    double temp_d = ui->horizontalSlider_2->value() / 100.0;
    ui->label_a->setText(QString::number(temp_d));
    a = temp_d;
    Control(0);
}

void HH_Model_Qt::on_horizontalSlider_3_sliderMoved() {
    ui->horizontalSlider_3->setRange(-11, 110);

    double temp_d = ui->horizontalSlider_3->value() / 100.0;
    ui->label_b->setText(QString::number(temp_d));
    b = temp_d;

    Control(0);
}

void HH_Model_Qt::on_horizontalSlider_4_sliderMoved() {
    ui->horizontalSlider_4->setRange(-70, 10);

    auto temp = ui->horizontalSlider_4->value();
    ui->label_c->setText(QString::number(temp));
    c = temp;

    Control(0);
}

void HH_Model_Qt::on_horizontalSlider_5_sliderMoved() {
    ui->horizontalSlider_5->setRange(-23, 41);

    auto temp = ui->horizontalSlider_5->value();
    ui->label_d->setText(QString::number(temp));
    d = temp;

    Control(0);
}

void HH_Model_Qt::on_pushButton_clicked() {
    a = 0.02;
    ui->horizontalSlider_2->setValue(a * 100);
    ui->label_a->setText(QString::number(a));

    b = 0.2;
    ui->horizontalSlider_3->setValue(b * 100);
    ui->label_b->setText(QString::number(b));

    c = -65.0;
    ui->horizontalSlider_4->setValue(c);
    ui->label_c->setText(QString::number(c));

    d = 6.0;
    ui->horizontalSlider_5->setValue(d);
    ui->label_d->setText(QString::number(d));

    num = T / 8 / dt;
    ui->horizontalSlider->setValue(T / 8);
    ui->label_time->setText(QString::number( T / 8));

    Control(0);
}
