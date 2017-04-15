#ifndef IZHIKEVICH_SPNET_H
#define IZHIKEVICH_SPNET_H

#include <QMainWindow>
#include "qcustomplot.h"
#include <utility>
#include <vector>

using std::vector;

namespace Ui {
class izhikevich_SPNET;
}

class izhikevich_SPNET : public QMainWindow
{
    Q_OBJECT

public:
    explicit izhikevich_SPNET(QWidget *parent = 0);
    ~izhikevich_SPNET();

    void Control(std::pair<QVector<double>, QVector<double>> temp_Pair);
    void PrintPlot(QCustomPlot *SPNET_Plot, std::pair<QVector<double>, QVector<double>> temp_Pair);
    long long Get_Random(const int &max);
    void Initialize();
    void Save_Data();
    std::pair<QVector<double>, QVector<double>> Data();
    void Simulation();

private slots:
    void on_pushButton_clicked();

private:
    Ui::izhikevich_SPNET *ui;

    const int Ne = 800; //excitatory neurons
    const int Ni = 200; // inhibitory neurons
    const int N = Ne + Ni; // total number of neurons
    const int M = 100; // the number of synapses per neuron
    const int D = 20; // maximal axonal conduction delay
    double sm = 10.0; // maximal synaptic strength
    vector<vector<long long>> post; //[N][M] indeces of postsynaptic neurons
    vector<vector<double>> s, sd; //[N][M] matrix of synaptic weights and their derivatives
    vector<vector<int>> delays_length; //[N][D] distribution of delays
    vector<vector<vector<int>>> delays; //[N][D][M] arrangement of delays
    vector<int> N_pre; //[N] presynaptic information
    vector<vector<int>> I_pre, D_pre; //[N][3*M] presynaptic information
    vector<vector<double*>> s_pre, sd_pre; //[N][3*M] presynaptic weights
    vector<vector<double>> LTP; //[N][1001 + D] STDP functions
    vector<double> LTD; //[N] STDP functions
    vector<double> a, d, V, U; //[N] neuronal dynamics parameters and activity variables
    int N_firings; // the number of fired neurons
    const int N_firings_max = 100 * N; // upper limit on the number of fired neurons per sec
    vector<vector<int>> firings; //[N_frings_max][2] // indeces and timings of spikes
    vector<double> I;
};

#endif // IZHIKEVICH_SPNET_H
