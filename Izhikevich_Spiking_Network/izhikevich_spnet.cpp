#include "izhikevich_spnet.h"
#include "ui_izhikevich_spnet.h"
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <random>
#include <ctime>

using namespace std;

izhikevich_SPNET::izhikevich_SPNET(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::izhikevich_SPNET)
{
    ui->setupUi(this);

    Initialize();

    ui->spinBox->setMinimum(1);
    ui->spinBox->setMaximum(864000);
}

izhikevich_SPNET::~izhikevich_SPNET()
{
    delete ui;
}

void izhikevich_SPNET::Control(pair<QVector<double>, QVector<double>> temp_Pair){
    ui->SPNET_Plot->clearGraphs();
    PrintPlot(ui->SPNET_Plot, temp_Pair);
    ui->SPNET_Plot->replot();
}

void izhikevich_SPNET::PrintPlot(QCustomPlot *SPNET_Plot, pair<QVector<double>, QVector<double>> temp_Pair){
    SPNET_Plot->addGraph(SPNET_Plot->xAxis, SPNET_Plot->yAxis);
    SPNET_Plot->graph(0)->addData(temp_Pair.first, temp_Pair.second);
    SPNET_Plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black), QBrush(Qt::black), 2));
    SPNET_Plot->graph(0)->setLineStyle(QCPGraph::lsNone);
    SPNET_Plot->graph(0)->rescaleAxes();
    SPNET_Plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}

void izhikevich_SPNET::Display_LCD(){
    ui->lcdNumber_2->setDecMode();
    ui->lcdNumber_2->setDigitCount(6);
    ui->lcdNumber_2->setSegmentStyle(QLCDNumber::Flat);
    ui->lcdNumber_2->display(count_sec);

    ui->lcdNumber->setDecMode();
    ui->lcdNumber->setDigitCount(6);
    ui->lcdNumber->setSegmentStyle(QLCDNumber::Flat);
    ui->lcdNumber->display(N_firings);
}

long long izhikevich_SPNET::Get_Random(const int &max) {
    static default_random_engine e(time(nullptr));
    uniform_int_distribution<long long> u(0, max - 1);

    return u(e);
}

void izhikevich_SPNET::Save_Data_post() {
    ofstream out("synapses.dat");
    for (auto i : post) {
        for(auto j : i){
            out << j << " ";
        }
        out << endl;
    }
}

void izhikevich_SPNET::Save_Data_firings() {
    ofstream out("spikes.dat");
    for (auto i = 1; i != N_firings; ++i) {
        if (firings[i][0] >= 0) {
            out << firings[i][0] << " " << firings[i][1] << endl;
        }
    }
}

pair<QVector<double>, QVector<double>> izhikevich_SPNET::Data_firings(){
    QVector<double> x, y0;

    for (auto i = 1; i != N_firings; ++i) {
        if (firings[i][0] >= 0) {
            x.push_back(firings[i][0]);
            y0.push_back(firings[i][1]);
        }
    }

    return make_pair(x, y0);
}

pair<QVector<double>, QVector<double>> izhikevich_SPNET::Data_N_firings(){
    vec_Nf_x.push_back(sec);
    vec_Nf_y0.push_back(N_firings);

    return make_pair(vec_Nf_x, vec_Nf_y0);
}

void izhikevich_SPNET::Initialize() {
    post.resize(N);
    for (auto &i : post) { i.resize(M); }
    s.resize(N);
    for (auto &i : s) { i.resize(M); }
    sd.resize(N);
    for (auto &i : sd) { i.resize(M); }
    delays_length.resize(N);
    for (auto &i : delays_length) { i.resize(D); }
    delays.resize(N);
    for (auto &i : delays) { i.resize(D); }
    for (auto &i : delays) { for (auto &j : i) { j.resize(M); } }
    N_pre.resize(N);
    I_pre.resize(N);
    for (auto &i : I_pre) { i.resize(3 * M); }
    D_pre.resize(N);
    for (auto &i : D_pre) { i.resize(3 * M); }
    s_pre.resize(N);
    for (auto &i : s_pre) { i.resize(3 * M); }
    sd_pre.resize(N);
    for (auto &i : sd_pre) { i.resize(3 * M); }
    LTP.resize(N);
    for (auto &i : LTP) { i.resize(1001 + D); }
    LTD.resize(N);
    a.resize(N);
    d.resize(N);
    V.resize(N);
    U.resize(N);
    firings.resize(N_firings_max);
    for (auto &i : firings) { i.resize(2); }
    I.resize(N);

    for (auto i = 0; i != Ne; ++i) a[i] = 0.02;// RS type
    for (auto i = Ne; i !=N; ++i) a[i] = 0.1; // FS type

    for (auto i = 0; i !=Ne; ++i) d[i] = 8.0; // RS type
    for (auto i = Ne; i != N; ++i) d[i] = 2.0; // FS type

    long long r;
    int exists;
    for (auto i = 0; i != N; ++i) {
        for (auto j = 0; j != M; ++j) {
            do {
                exists = 0; // avoid multiple synapses
                if (i < Ne) {
                    r = Get_Random(N);
                }
                else {
                    r = Get_Random(Ne); // inh -> exc only
                }
                if (r == i) {
                    exists = 1; // no self-synapses
                }
                for (auto k = 0; k != j; ++k) {
                    if (post[i][k] == r) {
                        exists = 1; // synapse already exists
                    }
                }
            } while (exists == 1);
            post[i][j] = r;
        }
    }

    for (auto i = 0; i != Ne; ++i)	for (auto j = 0; j != M; ++j) s[i][j] = 6.0;  // initial exc. synaptic weights
    for (auto i = Ne; i != N; ++i)	for (auto j = 0; j != M; ++j) s[i][j] = -5.0; // inhibitory synaptic weights
    for (auto i = 0; i != N; ++i)	for (auto j = 0; j != M; ++j) sd[i][j] = 0.0; // synaptic derivatives
    for (auto i = 0; i != N; ++i) {
        auto ind = 0;
        if (i < Ne) {
            for (auto j = 0; j != D; ++j) {
                delays_length[i][j] = M / D;	// uniform distribution of exc. synaptic delays
                for (auto k = 0; k != delays_length[i][j]; ++k) delays[i][j][k] = ind++;
            }
        }
        else {
            for (auto j = 0; j != D; ++j) delays_length[i][j] = 0;
            delays_length[i][0] = M;			// all inhibitory delays are 1 ms
            for (auto k = 0; k != delays_length[i][0]; ++k) delays[i][0][k] = ind++;
        }
    }

    for (auto i = 0; i != N; ++i) {
        N_pre[i] = 0;
        for (auto j = 0; j != Ne; ++j) {
            for (auto k = 0; k != M; ++k) {
                if (post[j][k] == i) {		// find all presynaptic neurons
                    I_pre[i][N_pre[i]] = j;	// add this neuron to the list
                    for (auto dd = 0; dd != D; ++dd){	// find the delay
                        for (auto jj = 0; jj != delays_length[j][dd]; ++jj){
                            if (post[j][delays[j][dd][jj]] == i) D_pre[i][N_pre[i]] = dd;
                        }
                    }
                    s_pre[i][N_pre[i]] = &s[j][k];	// pointer to the synaptic weight
                    sd_pre[i][N_pre[i]++] = &sd[j][k];// pointer to the derivative
                }
            }
        }
    }

    for (auto i = 0; i != N; ++i)   for (auto j = 0; j != 1 + D; ++j) LTP[i][j] = 0.0;
    for (auto i = 0; i != N; ++i)	LTD[i] = 0.0;
    for (auto i = 0; i != N; ++i)	V[i] = -65.0;		// initial values for V
    for (auto i = 0; i != N; ++i)	U[i] = 0.2 * V[i];	// initial values for U

    N_firings = 1;		// spike timings
    firings[0][0] = -D;	// put a dummy spike at -D for simulation efficiency
    firings[0][1] = 0;	// index of the dummy spike

    Save_Data_post();
}

void izhikevich_SPNET::Simulation() {
    for (; sec != T; sec += 1) {	// simulation of T sec
        for (auto t = 0; t != 1000; ++t) {			// simulation of 1 sec
            for (auto i = 0; i != N; ++i) I[i] = 0.0;	// reset the input
            for (auto i = 0; i != N / 1000; ++i) I[Get_Random(N)] = 20.0;		// random thalamic input

            for (auto i = 0; i != N; ++i) {
                if (V[i] >= 30) {					// did it fire?
                    V[i] = -65.0;					// voltage reset
                    U[i] += d[i];					// recovery variable reset
                    LTP[i][t + D] = 0.1;
                    LTD[i] = 0.12;
                    for (auto j = 0; j != N_pre[i]; ++j) *sd_pre[i][j] += LTP[I_pre[i][j]][t + D - D_pre[i][j] - 1];// this spike was after pre-synaptic spikes
                    firings[N_firings][0] = t;
                    firings[N_firings++][1] = i;
                    if (N_firings == N_firings_max) { cout << "Two many spikes at t=" << t << " (ignoring all)"; N_firings = 1; }
                }
            }

            auto k = N_firings;
            while (t - firings[--k][0] < D) {
                for (auto j = 0; j != delays_length[firings[k][1]][t - firings[k][0]]; ++j) {
                    auto i = post[firings[k][1]][delays[firings[k][1]][t - firings[k][0]][j]];
                    I[i] += s[firings[k][1]][delays[firings[k][1]][t - firings[k][0]][j]];
                    if (firings[k][1] < Ne) { // this spike is before postsynaptic spikes
                        sd[firings[k][1]][delays[firings[k][1]][t - firings[k][0]][j]] -= LTD[i];
                    }
                }
            }

            for (auto i = 0; i != N; ++i) {
                V[i] += 0.5*((0.04*V[i] + 5)*V[i] + 140 - U[i] + I[i]); // for numerical stability
                V[i] += 0.5*((0.04*V[i] + 5)*V[i] + 140 - U[i] + I[i]); // time step is 0.5 ms
                U[i] += a[i] * (0.2*V[i] - U[i]);
                LTP[i][t + D + 1] = 0.95*LTP[i][t + D];
                LTD[i] *= 0.95;
            }
        }

        ++count_sec; //Sum of sec;

        //cout << "sec=" << sec << ", firing rate=" << static_cast<double>(N_firings) / N << endl;

        //Save_Data_firings();

        Control(Data_firings());

        //Control(Data_N_firings());

        Display_LCD();

        for (auto i = 0; i != N; ++i) {		// prepare for the next sec
            for (auto j = 0; j != D + 1; ++j)
                LTP[i][j] = LTP[i][1000 + j];
        }
        auto k = N_firings - 1;
        while (1000 - firings[k][0] < D) k--;
        for (auto i = 1; i != N_firings - k; ++i) {
            firings[i][0] = firings[k + i][0] - 1000;
            firings[i][1] = firings[k + i][1];
        }
        N_firings -= k;

        for (auto i = 0; i != Ne; ++i) {	// modify only exc connections
            for (auto j = 0; j != M; ++j) {
                s[i][j] += 0.01 + sd[i][j];
                sd[i][j] *= 0.9;
                if (s[i][j] > sm) s[i][j] = sm;
                if (s[i][j] < 0) s[i][j] = 0.0;
            }
        }
    }
}

void izhikevich_SPNET::on_pushButton_clicked() {
    if(T < count_sec){
        vec_Nf_x.clear();
        vec_Nf_y0.clear();

        count_sec = 0.0;
        sec = 0;
        Initialize();
        Simulation();
    }
    else{
        Simulation();
    }
}

void izhikevich_SPNET::on_spinBox_editingFinished() {
    T = ui->spinBox->value();
}

void izhikevich_SPNET::on_pushButton_2_clicked() {
    ++T;
    on_pushButton_clicked();
    ui->spinBox->setValue(T);
    T = ui->spinBox->value();
}
