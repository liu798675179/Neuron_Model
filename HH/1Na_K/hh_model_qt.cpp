#include "hh_model_qt.h"
#include "ui_hh_model_qt.h"
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
        ui->HH_Plot1->clearGraphs();
        PrintPlot(ui->HH_Plot, ui->HH_Plot1);
        setWindowTitle(WinName);
        ui->HH_Plot->replot();
        ui->HH_Plot1->replot();
    }
    else if(cmd == 1){
        ui->HH_Plot->clearGraphs();
        ui->HH_Plot1->clearGraphs();
        PrintPlot1(ui->HH_Plot, ui->HH_Plot1);
        setWindowTitle(WinName);
        ui->HH_Plot->replot();
        ui->HH_Plot1->replot();
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

auto num = 5000;

void HH_Model_Qt::PrintPlot(QCustomPlot *HH_Plot, QCustomPlot *HH_Plot1){
    WinName = "The INa + IKdr-minimal model.";

    auto t = 0.0;
    QVector<double> x(num), y0(num), y1(num), y2(num);

    Neuron_Normal temp_Neuron_Normal;
    auto Tuple_Normal = temp_Neuron_Normal.Simulation();

    Neuron_AD temp_Neuron_AD;
    auto Tuple_AD = temp_Neuron_AD.Simulation();

    // add two new graphs and set their look:
    HH_Plot->addGraph(HH_Plot->xAxis, HH_Plot->yAxis);
    HH_Plot->graph(0)->setPen(QPen(Qt::gray));
    HH_Plot->graph(0)->setName("The membrane potential.");
    HH_Plot->addGraph(HH_Plot->xAxis, HH_Plot->yAxis2);
    HH_Plot->graph(1)->setPen(QPen(Qt::red));
    HH_Plot->graph(1)->setName("The current of Na.");
    HH_Plot->addGraph(HH_Plot->xAxis, HH_Plot->yAxis2);
    HH_Plot->graph(2)->setPen(QPen(Qt::green));
    HH_Plot->graph(2)->setName("The current of K.");
    HH_Plot->legend->setVisible(true);

    // generate some points of data (y0 for first, y1 for second graph):
    for(auto i = 0; i < num; ++i){
        x[i] = t;
        y0[i] = get<0>(Tuple_Normal)[i];
        y1[i] = get<1>(Tuple_Normal)[i];
        y2[i] = get<2>(Tuple_Normal)[i];

        t += dt;
    }

    // configure right and top axis to show ticks but no labels:
    HH_Plot->xAxis->setLabel("Time: mS");
    HH_Plot->xAxis2->setLabel("Normal");
    HH_Plot->yAxis->setLabel("Potential: mV");
    HH_Plot->yAxis2->setLabel("Current: nA");
    HH_Plot->xAxis2->setVisible(true);
    HH_Plot->xAxis2->setTickLabels(false);
    HH_Plot->yAxis2->setVisible(true);
    HH_Plot->yAxis2->setTickLabels(true);

    // make left and bottom axes always transfer their ranges to right and top axes:
    connect(HH_Plot->xAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot->xAxis2, SLOT(setRange(QCPRange)));
    connect(HH_Plot->yAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot->yAxis2, SLOT(setRange(QCPRange)));

    // pass data points to graphs:
    HH_Plot->graph(0)->addData(x, y0);
    HH_Plot->graph(1)->addData(x, y1);
    HH_Plot->graph(2)->addData(x, y2);

    // let the ranges scale themselves so graph 0 fits perfectly in the visible area:
    HH_Plot->graph(0)->rescaleAxes();
    // same thing for graph 1, but only enlarge ranges (in case graph 1 is smaller than graph 0):
    HH_Plot->graph(1)->rescaleAxes(true);
    HH_Plot->graph(2)->rescaleAxes(true);

    // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
    HH_Plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);

    HH_Plot1->addGraph(HH_Plot1->xAxis, HH_Plot1->yAxis);
    HH_Plot1->graph(0)->setPen(QPen(Qt::gray));
    HH_Plot1->graph(0)->setName("The membrane potential.");
    HH_Plot1->addGraph(HH_Plot1->xAxis, HH_Plot1->yAxis2);
    HH_Plot1->graph(1)->setPen(QPen(Qt::red));
    HH_Plot1->graph(1)->setName("The current of Na.");
    HH_Plot1->addGraph(HH_Plot1->xAxis, HH_Plot1->yAxis2);
    HH_Plot1->graph(2)->setPen(QPen(Qt::green));
    HH_Plot1->graph(2)->setName("The current of K.");
    HH_Plot1->legend->setVisible(true);

    //reset some variables
    x.clear();
    x.resize(num);
    y0.clear();
    y0.resize(num);
    y1.clear();
    y1.resize(num);
    t = 0;

    for (auto i = 0; i < num; ++i) {
        x[i] = t;
        y0[i] = get<0>(Tuple_AD)[i];
        y1[i] = get<1>(Tuple_AD)[i];
        y2[i] = get<2>(Tuple_AD)[i];

        t += dt;
    }

    HH_Plot1->xAxis->setLabel("Time: mS");
    HH_Plot1->xAxis2->setLabel("AD");
    HH_Plot1->yAxis->setLabel("Potential: mV");
    HH_Plot1->yAxis2->setLabel("Current: nA");
    HH_Plot1->xAxis2->setVisible(true);
    HH_Plot1->xAxis2->setTickLabels(false);
    HH_Plot1->yAxis2->setVisible(true);
    HH_Plot1->yAxis2->setTickLabels(true);

    connect(HH_Plot1->xAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot1->xAxis2, SLOT(setRange(QCPRange)));
    connect(HH_Plot1->yAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot1->yAxis2, SLOT(setRange(QCPRange)));

    HH_Plot1->graph(0)->addData(x, y0);
    HH_Plot1->graph(1)->addData(x, y1);
    HH_Plot1->graph(2)->addData(x, y2);

    HH_Plot1->graph(0)->rescaleAxes();
    HH_Plot1->graph(1)->rescaleAxes(true);
    HH_Plot1->graph(2)->rescaleAxes(true);

    HH_Plot1->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
}

void HH_Model_Qt::PrintPlot1(QCustomPlot *HH_Plot, QCustomPlot *HH_Plot1){
    WinName = "The INa + IKdr-minimal model.";

    auto t = 0.0;
    QVector<double> x(num), y0(num), y1(num), y2(num);

    Neuron_Normal temp_Neuron_Normal;
    auto Tuple_Normal = temp_Neuron_Normal.Simulation();

    Neuron_AD temp_Neuron_AD;
    auto Tuple_AD = temp_Neuron_AD.Simulation();

    // add two new graphs and set their look:
    HH_Plot->addGraph(HH_Plot->xAxis, HH_Plot->yAxis);
    HH_Plot->graph(0)->setPen(QPen(Qt::gray));
    HH_Plot->graph(0)->setName("The membrane potential.");
    HH_Plot->addGraph(HH_Plot->xAxis, HH_Plot->yAxis2);
    HH_Plot->graph(1)->setPen(QPen(Qt::red));
    HH_Plot->graph(1)->setName("The current of -Na.");
    HH_Plot->addGraph(HH_Plot->xAxis, HH_Plot->yAxis2);
    HH_Plot->graph(2)->setPen(QPen(Qt::green));
    HH_Plot->graph(2)->setName("The current of K.");
    HH_Plot->legend->setVisible(true);

    // generate some points of data (y0 for first, y1 for second graph):
    for(auto i = 0; i < num; ++i){
        x[i] = t;
        y0[i] = get<0>(Tuple_Normal)[i];
        y1[i] = -get<1>(Tuple_Normal)[i];
        y2[i] = get<2>(Tuple_Normal)[i];

        t += dt;
    }

    // configure right and top axis to show ticks but no labels:
    HH_Plot->xAxis->setLabel("Time: mS");
    HH_Plot->xAxis2->setLabel("Normal");
    HH_Plot->yAxis->setLabel("Potential: mV");
    HH_Plot->yAxis2->setLabel("Current: nA");
    HH_Plot->xAxis2->setVisible(true);
    HH_Plot->xAxis2->setTickLabels(false);
    HH_Plot->yAxis2->setVisible(true);
    HH_Plot->yAxis2->setTickLabels(true);

    // make left and bottom axes always transfer their ranges to right and top axes:
    connect(HH_Plot->xAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot->xAxis2, SLOT(setRange(QCPRange)));
    connect(HH_Plot->yAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot->yAxis2, SLOT(setRange(QCPRange)));

    // pass data points to graphs:
    HH_Plot->graph(0)->addData(x, y0);
    HH_Plot->graph(1)->addData(x, y1);
    HH_Plot->graph(2)->addData(x, y2);

    // let the ranges scale themselves so graph 0 fits perfectly in the visible area:
    HH_Plot->graph(0)->rescaleAxes();
    // same thing for graph 1, but only enlarge ranges (in case graph 1 is smaller than graph 0):
    HH_Plot->graph(1)->rescaleAxes(true);
    HH_Plot->graph(2)->rescaleAxes(true);

    // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
    HH_Plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);

    HH_Plot1->addGraph(HH_Plot1->xAxis, HH_Plot1->yAxis);
    HH_Plot1->graph(0)->setPen(QPen(Qt::gray));
    HH_Plot1->graph(0)->setName("The membrane potential.");
    HH_Plot1->addGraph(HH_Plot1->xAxis, HH_Plot1->yAxis2);
    HH_Plot1->graph(1)->setPen(QPen(Qt::red));
    HH_Plot1->graph(1)->setName("The current of -Na.");
    HH_Plot1->addGraph(HH_Plot1->xAxis, HH_Plot1->yAxis2);
    HH_Plot1->graph(2)->setPen(QPen(Qt::green));
    HH_Plot1->graph(2)->setName("The current of K.");
    HH_Plot1->legend->setVisible(true);

    //reset some variables
    x.clear();
    x.resize(num);
    y0.clear();
    y0.resize(num);
    y1.clear();
    y1.resize(num);
    t = 0;

    for (auto i = 0; i < num; ++i) {
        x[i] = t;
        y0[i] = get<0>(Tuple_AD)[i];
        y1[i] = -get<1>(Tuple_AD)[i];
        y2[i] = get<2>(Tuple_AD)[i];

        t += dt;
    }

    HH_Plot1->xAxis->setLabel("Time: mS");
    HH_Plot1->xAxis2->setLabel("AD");
    HH_Plot1->yAxis->setLabel("Potential: mV");
    HH_Plot1->yAxis2->setLabel("Current: nA");
    HH_Plot1->xAxis2->setVisible(true);
    HH_Plot1->xAxis2->setTickLabels(false);
    HH_Plot1->yAxis2->setVisible(true);
    HH_Plot1->yAxis2->setTickLabels(true);

    connect(HH_Plot1->xAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot1->xAxis2, SLOT(setRange(QCPRange)));
    connect(HH_Plot1->yAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot1->yAxis2, SLOT(setRange(QCPRange)));

    HH_Plot1->graph(0)->addData(x, y0);
    HH_Plot1->graph(1)->addData(x, y1);
    HH_Plot1->graph(2)->addData(x, y2);

    HH_Plot1->graph(0)->rescaleAxes();
    HH_Plot1->graph(1)->rescaleAxes(true);
    HH_Plot1->graph(2)->rescaleAxes(true);

    HH_Plot1->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
}

void HH_Model_Qt::on_doubleSpinBox_editingFinished() {
    if(ui->doubleSpinBox->value() != 0) {
        Iapp = ui->doubleSpinBox->text().toDouble();
    }
}

void HH_Model_Qt::on_horizontalSlider_sliderPressed() {
    ui->label_3->setVisible(true);
    ui->label_3->setText(QString::number(ui->horizontalSlider->value()) + " mS");
    num = ui->horizontalSlider->value() / dt;
}

void HH_Model_Qt::on_horizontalSlider_sliderMoved() {
    ui->horizontalSlider->setRange(5, 1000);
    ui->label_3->setVisible(true);
    ui->label_3->setText(QString::number(ui->horizontalSlider->value()) + " mS");
    num = ui->horizontalSlider->value() / dt;
}

void HH_Model_Qt::on_pushButton_clicked() {
    Control(0);
}

int count = 1;
void HH_Model_Qt::on_pushButton_2_clicked() {
    if(count % 2 == 0){
        Control(0);
        ++count;
    }else{
        Control(1);
        ++count;
    }
}
