#ifndef IZHIKEVICH_SPNET_H
#define IZHIKEVICH_SPNET_H

#include <QMainWindow>
#include "qcustomplot.h"
#include <utility>
#include <vector>
#include <map>

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

    void Control(int it, std::pair<QVector<double>, QVector<double>> temp_Pair);
    void PrintPlot(QCustomPlot *SPNET_Plot, std::pair<QVector<double>, QVector<double>> temp_Pair);
    void PrintPlot_2(QCustomPlot *SPNET_Plot_2, std::pair<QVector<double>, QVector<double>> temp_Pair);
    void PrintPlot_3(QCustomPlot *SPNET_Plot_3, std::pair<QVector<double>, QVector<double>> temp_Pair);
    int Get_Random_int(const int &max);
    int Get_Random_int(const int &min, const int &max);
    double Get_Random_real();
    double Noise();
    void Initialize();
    void Save_Data_post();
    void Save_Data_delay();
    void Save_Data_weight();
    void Save_Data_current();
    void Save_Data_firings();
    void Display_LCD();
    std::pair<QVector<double>, QVector<double>> Data_firings();
    std::pair<QVector<double>, QVector<double>> Data_N_firings();
	std::pair<QVector<double>, QVector<double>> Data_Power_Law();
    void Simulation();

private slots:
    void on_pushButton_clicked();

    void on_spinBox_editingFinished();

    void on_pushButton_2_clicked();

private:
    Ui::izhikevich_SPNET *ui;

    const int Ne = 800; //excitatory neurons
    const int Ni = 200; // inhibitory neurons
    const int N = Ne + Ni; // total number of neurons
    const int M = 100; // the number of synapses per neuron
    const int D = 20; // maximal axonal conduction delay
    double sm = 10.0; // maximal synaptic strength
    vector<vector<int>> post; //[N][M] indeces of postsynaptic neurons
    vector<vector<double>> s, sd; //[N][M] matrix of synaptic weights and their derivatives
    vector<vector<int>> delays_length; //[N][D] distribution of delays
    vector<vector<vector<int>>> delays; //[N][D][M] arrangement of delays
    vector<int> N_pre; //[N] presynaptic information
    vector<vector<int>> I_pre, D_pre; //[N][3*M] presynaptic information
    vector<vector<double*>> sd_pre; //[N][3*M] presynaptic weights
    vector<vector<double>> LTP; //[N][1001 + D] STDP functions
    vector<double> LTD; //[N] STDP functions
    vector<double> a, b, c, d, V, U; //[N] neuronal dynamics parameters and activity variables
    int N_firings; // the number of fired neurons
    const int N_firings_max = 100 * N; // upper limit on the number of fired neurons per sec
    vector<vector<int>> firings; //[N_frings_max][2] // indeces and timings of spikes
    vector<double> I; //Current
    vector<vector<double>> Save_I; // Current of 1000ms
    long long sec = 0; //Simulation time
    long long T = 1; //Total simulaton time
    double count_sec = 0.0; //Sum of sec
    int temp_10 = 0;
    QVector<double> vec_f_x, vec_f_y0, vec_f_x_10; //Save sec and firings
    QVector<double> vec_Nf_x, vec_Nf_y0; //Save sec and N_firings
	std::map<int, int> map_size_count;
};

#endif // IZHIKEVICH_SPNET_H
