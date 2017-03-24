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
        ui->HH_Plot->clearPlottables();
        ui->HH_Plot1->clearPlottables();
        PrintPlot(ui->HH_Plot, ui->HH_Plot1);
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

auto pointcount = T / dt;

void HH_Model_Qt::PrintPlot(QCustomPlot *HH_Plot, QCustomPlot *HH_Plot1){
    WinName = "The INa + IKdr-minimal model.";

    QVector<double> x, x1, y0, y1;
    QVector<QCPCurveData> phaseorbit;

    Neuron_Normal temp_Neuron_Normal;
    auto n_null_Normal = temp_Neuron_Normal.n_nullcline();
    auto V_null_Normal = temp_Neuron_Normal.V_nullcline();
    auto P_O_Normal = temp_Neuron_Normal.Phase_Orbit();

    Neuron_AD temp_Neuron_AD;
    auto n_null_AD = temp_Neuron_AD.n_nullcline();
    auto V_null_AD = temp_Neuron_AD.V_nullcline();
    auto P_O_AD = temp_Neuron_AD.Phase_Orbit();

    // add two new graphs and set their look:
    HH_Plot->addGraph(HH_Plot->xAxis, HH_Plot->yAxis);
    HH_Plot->graph(0)->setPen(QPen(Qt::red));
    HH_Plot->graph(0)->setName("n-nullcline");
    HH_Plot->addGraph(HH_Plot->xAxis, HH_Plot->yAxis);
    HH_Plot->graph(1)->setPen(QPen(Qt::green));
    HH_Plot->graph(1)->setName("V-nullcline");
    QCPCurve *PhaseOrbit = new QCPCurve(HH_Plot->xAxis, HH_Plot->yAxis);
    PhaseOrbit->setPen(QPen(Qt::black));
    PhaseOrbit->setName("Phase Orbit");
    HH_Plot->legend->setVisible(true);
    HH_Plot->legend->setBrush(QColor(0,0,0,0));

    // generate some points of data (y0 for first, y1 for second graph):
    for(auto i = -80.0; i < 60.0; i += dt){
        x.push_back(i);
        x1.push_back(i);
    }
    for(auto &i : n_null_Normal){
        y0.push_back(i);
    }
    for(auto &i : V_null_Normal){
        y1.push_back(i);
    }

    for(auto i = 0; i <= pointcount; ++i){
        phaseorbit.push_back(QCPCurveData(i, P_O_Normal.first[i], P_O_Normal.second[i]));
    }

    // configure right and top axis to show ticks but no labels:
    HH_Plot->xAxis->setLabel("V/mV");
    HH_Plot->xAxis2->setLabel("Normal");
    HH_Plot->yAxis->setLabel("n");
    HH_Plot->xAxis2->setVisible(true);
    HH_Plot->xAxis2->setTickLabels(false);
    HH_Plot->yAxis2->setVisible(true);
    HH_Plot->yAxis2->setTickLabels(false);

    // make left and bottom axes always transfer their ranges to right and top axes:
    connect(HH_Plot->xAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot->xAxis2, SLOT(setRange(QCPRange)));
    connect(HH_Plot->yAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot->yAxis2, SLOT(setRange(QCPRange)));

    // pass data points to graphs:
    HH_Plot->graph(0)->addData(x, y0);
    HH_Plot->graph(1)->addData(x1, y1);
    PhaseOrbit->data()->set(phaseorbit, true);

    // let the ranges scale themselves so graph 0 fits perfectly in the visible area:
    HH_Plot->graph(0)->rescaleAxes();
    // same thing for graph 1, but only enlarge ranges (in case graph 1 is smaller than graph 0):
    HH_Plot->graph(1)->rescaleAxes(true);

    // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
    HH_Plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    HH_Plot1->addGraph(HH_Plot1->xAxis, HH_Plot1->yAxis);
    HH_Plot1->graph(0)->setPen(QPen(Qt::red));
    HH_Plot1->graph(0)->setName("n-nullcline");
    HH_Plot1->addGraph(HH_Plot1->xAxis, HH_Plot1->yAxis);
    HH_Plot1->graph(1)->setPen(QPen(Qt::green));
    HH_Plot1->graph(1)->setName("V-nullcline");
    QCPCurve *PhaseOrbit1 = new QCPCurve(HH_Plot1->xAxis, HH_Plot1->yAxis);
    PhaseOrbit1->setPen(QPen(Qt::black));
    PhaseOrbit1->setName("Phase Orbit");
    HH_Plot1->legend->setVisible(true);
    HH_Plot1->legend->setBrush(QColor(0,0,0,0));

    //reset some variables
    y0.clear();
    y1.clear();
    phaseorbit.clear();

    for(auto &i : n_null_AD){
        y0.push_back(i);
    }
    for(auto &i : V_null_AD){
        y1.push_back(i);
    }

    for(auto i = 0; i <= pointcount; ++i){
        phaseorbit.push_back(QCPCurveData(i, P_O_AD.first[i], P_O_AD.second[i]));
    }

    HH_Plot1->xAxis->setLabel("V/mV");
    HH_Plot1->xAxis2->setLabel("AD");
    HH_Plot1->yAxis->setLabel("n");
    HH_Plot1->xAxis2->setVisible(true);
    HH_Plot1->xAxis2->setTickLabels(false);
    HH_Plot1->yAxis2->setVisible(true);
    HH_Plot1->yAxis2->setTickLabels(false);

    connect(HH_Plot1->xAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot1->xAxis2, SLOT(setRange(QCPRange)));
    connect(HH_Plot1->yAxis, SIGNAL(rangeChanged(QCPRange)), HH_Plot1->yAxis2, SLOT(setRange(QCPRange)));

    HH_Plot1->graph(0)->addData(x, y0);
    HH_Plot1->graph(1)->addData(x1, y1);
    PhaseOrbit1->data()->set(phaseorbit, true);

    HH_Plot1->graph(0)->rescaleAxes();
    HH_Plot1->graph(1)->rescaleAxes(true);

    HH_Plot1->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

void HH_Model_Qt::on_doubleSpinBox_editingFinished() {
    if(ui->doubleSpinBox->value() >= 0) {
        Iapp = ui->doubleSpinBox->text().toDouble();
    }
}

void HH_Model_Qt::on_pushButton_clicked() {
    Control(0);
}
