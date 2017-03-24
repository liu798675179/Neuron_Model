#ifndef _NEURON_H_
#define _NEURON_H_

#include <cmath>
#include <vector>
#include <utility>
#include <iostream>
#include <tuple>

using std::exp;
using std::pow;
using std::get;
using std::cout;
using std::endl;
using std::pair;
using std::tuple;
using std::vector;

/*************************************************************
**********************Some constants value.*******************
*************************************************************/

//The membrane capacitance. Demention : ¦ÌF/cm^2
double const C = 1;
//The current conductance of k+. Demention : mS/cm^2
double g_Kdr;
//The current conductance of Na+. Demention : mS/cm^2
double const g_Na = 35;
//The current conductance of other ions. Demention : mS/cm^2
double const g_L = 0.05;
//The membrane potential of K+. Demention : mV
double const V_k = -90;
//The membrane potential of Na+. Demention : mV
double const V_Na = 55;
//The membrane potential of other ions. Demention : mV
double const V_L = -70;
//The step length.
double const dt = 0.001;
//The stimulation current. Demention : nA
double Iapp = 6;
//The simulation times. Demention : mS
double T = 1000.0;

//The INa + IKdr-minimal model.
class Neuron_Normal {
public:
    Neuron_Normal() = default;
    Neuron_Normal(double const &temp_v, double const &temp_h, double const &temp_n) : V(temp_v), h(temp_h), n(temp_n) { }
    Neuron_Normal(double const &temp_v, double const &temp_I_Na, double const &temp_I_K, double const &temp_h, double const &temp_n) : V(temp_v), I_Na(temp_I_Na), I_K(temp_I_K), h(temp_h), n(temp_n) { }
    virtual ~Neuron_Normal() = default;

    void initialize(double const &, double const &, double const &);

    double n_V(Neuron_Normal &) const;
    double h_v(Neuron_Normal &) const;
    double m_V(Neuron_Normal &) const;
    double tau_hV(Neuron_Normal &) const;
    double tau_nV(Neuron_Normal &) const;
    Neuron_Normal HH_Model(Neuron_Normal &) const;

    double Energy(Neuron_Normal &) const;
    tuple<vector<double>, vector<double>, vector<double>> Simulation() const;

protected:
    //The membrane potential.
    double V = -65;
    //The current of Na.
    double I_Na = 0;
    //The current of K.
    double I_K = 0;

    //The transient Na+ current inactivation variable h.
    double h = 0.1;
    //The delayed rectified K+ current activation variables n.
    double n = 0.1;
    //The energy consumption.
    double H = 0.0;
};

inline
void Neuron_Normal::initialize(double const &temp_v, double const &temp_h, double const &temp_n) {
    V = temp_v;
    h = temp_h;
    n = temp_n;
}

inline
double Neuron_Normal::n_V(Neuron_Normal &temp_Neu) const {
    return 1 / (1 + exp((-35 - temp_Neu.V) / 10));
}

inline
double Neuron_Normal::h_v(Neuron_Normal &temp_Neu) const {
    return 1 / (1 + exp((-45 - temp_Neu.V) / -7));
}

inline
double Neuron_Normal::m_V(Neuron_Normal &temp_Neu) const {
    return 1 / (1 + exp((-30 - temp_Neu.V) / 9.5));
}

inline
double Neuron_Normal::tau_hV(Neuron_Normal &temp_Neu) const {
    return 0.1 + 0.75 * 1 / (1 + exp((-40.5 - temp_Neu.V) / -6));
}

inline
double Neuron_Normal::tau_nV(Neuron_Normal &temp_Neu) const {
    return 0.1 + 0.5 * 1 / (1 + exp((-27 - temp_Neu.V) / -15));
}

inline
Neuron_Normal Neuron_Normal::HH_Model(Neuron_Normal &temp_Neu) const {
    g_Kdr = 6;
    temp_Neu.I_Na = g_Na * pow(m_V(temp_Neu), 3) * temp_Neu.h * (temp_Neu.V - V_Na);
    temp_Neu.I_K = g_Kdr * pow(temp_Neu.n, 4) * (temp_Neu.V - V_k);
    temp_Neu.V += dt * (-g_L * (temp_Neu.V - V_L) - g_Na * pow(m_V(temp_Neu), 3) * temp_Neu.h * (temp_Neu.V - V_Na) - g_Kdr * pow(temp_Neu.n, 4) * (temp_Neu.V - V_k) + Iapp);
    temp_Neu.h += dt * (h_v(temp_Neu) - temp_Neu.h) / tau_hV(temp_Neu);
    temp_Neu.n += dt * (n_V(temp_Neu) - temp_Neu.n) / tau_nV(temp_Neu);

    return Neuron_Normal(temp_Neu.V, temp_Neu.I_Na, temp_Neu.I_K, temp_Neu.h, temp_Neu.n);
}

inline
double Neuron_Normal::Energy(Neuron_Normal &temp_Neu) const {
    g_Kdr = 6;
    temp_Neu.H += dt * (temp_Neu.V * Iapp - g_L * pow((temp_Neu.V - V_L), 2) - g_Na * pow(m_V(temp_Neu), 3) * temp_Neu.h * pow((temp_Neu.V - V_Na), 2) - g_Kdr * pow(temp_Neu.n, 4) * pow((temp_Neu.V - V_k), 2));

    return temp_Neu.H;
}


inline
tuple<vector<double>, vector<double>, vector<double>> Neuron_Normal::Simulation() const {
    vector<double> Vec_Potential;
    vector<double> Vec_Na_Potential;
    vector<double> Vec_K_Potential;
    //vector<double> Vec_Energy;
    //double temp_E;
    Neuron_Normal temp_Neu;

    for (auto i = 0.0; i <= T; i += dt) {
        //temp_E = Energy(temp_Neu);
        temp_Neu = HH_Model(temp_Neu);

        Vec_Potential.push_back(temp_Neu.V);
        Vec_Na_Potential.push_back(temp_Neu.I_Na);
        Vec_K_Potential.push_back(temp_Neu.I_K);
        //Vec_Energy.push_back(temp_E);
    }

    //return std::make_pair(Vec_Potential, Vec_Energy);
    return std::make_tuple(Vec_Potential, Vec_Na_Potential, Vec_K_Potential);
}

//The INa + IKdr-minimal model.
class Neuron_AD {
public:
    Neuron_AD() = default;
    Neuron_AD(double const &temp_v, double const &temp_h, double const &temp_n) : V(temp_v), h(temp_h), n(temp_n) { }
    Neuron_AD(double const &temp_v, double const &temp_I_Na, double const &temp_I_K, double const &temp_h, double const &temp_n) : V(temp_v), I_Na(temp_I_Na), I_K(temp_I_K), h(temp_h), n(temp_n) { }
    virtual ~Neuron_AD() = default;

    void initialize(double const &, double const &, double const &);

    double n_V(Neuron_AD &) const;
    double h_v(Neuron_AD &) const;
    double m_V(Neuron_AD &) const;
    double tau_hV(Neuron_AD &) const;
    double tau_nV(Neuron_AD &) const;
    Neuron_AD HH_Model(Neuron_AD &) const;

    double Energy(Neuron_AD &) const;
    tuple<vector<double>, vector<double>, vector<double>> Simulation() const;

protected:
    //The membrane potential.
    double V = -65;
    //The current of Na.
    double I_Na = 0;
    //The current of Na.
    double I_K = 0;
    //The transient Na+ current inactivation variable h.
    double h = 0.1;
    //The delayed rectified K+ current activation variables n.
    double n = 0.1;
    //The energy consumption.
    double H = 0.0;
};

inline
void Neuron_AD::initialize(double const &temp_v, double const &temp_h, double const &temp_n) {
    V = temp_v;
    h = temp_h;
    n = temp_n;
}

inline
double Neuron_AD::n_V(Neuron_AD &temp_Neu) const {
    return 1 / (1 + exp((-43.1 - temp_Neu.V) / 7.3));
}

inline
double Neuron_AD::h_v(Neuron_AD &temp_Neu) const {
    return 1 / (1 + exp((-45 - temp_Neu.V) / -7));
}

inline
double Neuron_AD::m_V(Neuron_AD &temp_Neu) const {
    return 1 / (1 + exp((-30 - temp_Neu.V) / 9.5));
}

inline
double Neuron_AD::tau_hV(Neuron_AD &temp_Neu) const {
    return 0.1 + 0.75 * 1 / (1 + exp((-40.5 - temp_Neu.V) / -6));
}

inline
double Neuron_AD::tau_nV(Neuron_AD &temp_Neu) const {
    return 0.1 + 0.5 * 1 / (1 + exp((-27 - temp_Neu.V) / -15));
}

inline
Neuron_AD Neuron_AD::HH_Model(Neuron_AD &temp_Neu) const {
    g_Kdr = 10;
    temp_Neu.I_Na = g_Na * pow(m_V(temp_Neu), 3) * temp_Neu.h * (temp_Neu.V - V_Na);
    temp_Neu.I_K = g_Kdr * pow(temp_Neu.n, 4) * (temp_Neu.V - V_k);
    temp_Neu.V += dt * (-g_L * (temp_Neu.V - V_L) - g_Na * pow(m_V(temp_Neu), 3) * temp_Neu.h * (temp_Neu.V - V_Na) - g_Kdr * pow(temp_Neu.n, 4) * (temp_Neu.V - V_k) + Iapp);
    temp_Neu.h += dt * (h_v(temp_Neu) - temp_Neu.h) / tau_hV(temp_Neu);
    temp_Neu.n += dt * (n_V(temp_Neu) - temp_Neu.n) / tau_nV(temp_Neu);

    return Neuron_AD(temp_Neu.V, temp_Neu.I_Na, temp_Neu.I_K, temp_Neu.h, temp_Neu.n);
}

inline
double Neuron_AD::Energy(Neuron_AD &temp_Neu) const {
    g_Kdr = 10;
    temp_Neu.H += dt * (temp_Neu.V * Iapp - g_L * pow((temp_Neu.V - V_L), 2) - g_Na * pow(m_V(temp_Neu), 3) * temp_Neu.h * pow((temp_Neu.V - V_Na), 2) - g_Kdr * pow(temp_Neu.n, 4) * pow((temp_Neu.V - V_k), 2));

    return temp_Neu.H;
}


inline
tuple<vector<double>, vector<double>, vector<double>> Neuron_AD::Simulation() const {
    vector<double> Vec_Potential;
    vector<double> Vec_Na_Potential;
    vector<double> Vec_K_Potential;
    //vector<double> Vec_Energy;
    //double temp_E;
    Neuron_AD temp_Neu;

    for (auto i = 0.0; i <= T; i += dt) {
        //temp_E = Energy(temp_Neu);
        temp_Neu = HH_Model(temp_Neu);

        Vec_Potential.push_back(temp_Neu.V);
        Vec_Na_Potential.push_back(temp_Neu.I_Na);
        Vec_K_Potential.push_back(temp_Neu.I_K);
        //Vec_Energy.push_back(temp_E);
    }

    //return std::make_pair(Vec_Potential, Vec_Energy);
    return std::make_tuple(Vec_Potential, Vec_Na_Potential, Vec_K_Potential);
}

#endif _NEURON_H_
