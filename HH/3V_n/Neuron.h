#ifndef _NEURON_H_
#define _NEURON_H_

#include <cmath>
#include <vector>
#include <utility>
#include <iostream>
#include <tuple>
#include <random>
#include <ctime>

#include "Quartic_Equation.h"

using std::exp;
using std::pow;
using std::get;
using std::cout;
using std::endl;
using std::pair;
using std::tuple;
using std::vector;
using std::make_pair;

std::default_random_engine e(time(nullptr));
std::uniform_real_distribution<double> u;

/*************************************************************
**********************Some constants value.*******************
*************************************************************/

double const Pi = 3.1415926;

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
double Iapp = 0;
//The simulation times. Demention : mS
double T = 1000.0;

//The INa + IKdr-minimal model.
class Neuron_Normal {
public:
    Neuron_Normal() = default;
    Neuron_Normal(double const &temp_v, double const &temp_n) : V(temp_v),  n(temp_n) { }
    virtual ~Neuron_Normal() = default;

    void initialize(double const &, double const &);

    double n_V(double const &) const;
    double m_V(double const &) const;
    double tau_nV(double const &) const;
    Neuron_Normal HH_Model(Neuron_Normal &) const;

    double A(double const &) const;
    static double B(double const &);
    static double C(double const &);
    Quartic_Equation Solve_n(double const &) const;

    vector<double> n_nullcline() const;
    vector<double> V_nullcline() const;
    pair<vector<double>, vector<double>> Phase_Orbit() const;

    double Noise() const{
        return std::sqrt(std::abs(-4.0 * dt * std::log(u(e)) * std::cos(2 * Pi * u(e))));
    }

protected:
    //The membrane potential.
    double V = -65;

    //The delayed rectified K+ current activation variables n.
    double n = 0.8;
    //The transient Na+ current inactivation variable h, linear relation with n.
    double h = 0.8820 - 1.0028 * n;
};

inline
void Neuron_Normal::initialize(double const &temp_v, double const &temp_n) {
    V = temp_v;
    n = temp_n;
}

inline
double Neuron_Normal::n_V(double const &temp_V) const {
    return 1 / (1 + exp((-35 - temp_V) / 10));
}

inline
double Neuron_Normal::m_V(double const &temp_V) const {
    return 1 / (1 + exp((-30 - temp_V) / 9.5));
}


inline
double Neuron_Normal::tau_nV(double const &temp_V) const {
    return 0.1 + 0.5 * 1 / (1 + exp((-27 - temp_V) / -15));
}

inline
Neuron_Normal Neuron_Normal::HH_Model(Neuron_Normal &temp_Neu) const {
    g_Kdr = 6;

    temp_Neu.V += dt * (-g_L * (temp_Neu.V - V_L) - g_Na * pow(m_V(temp_Neu.V), 3) * temp_Neu.h * (temp_Neu.V - V_Na) - g_Kdr * pow(temp_Neu.n, 4) * (temp_Neu.V - V_k) + Iapp);// + Noise();
    temp_Neu.n += dt * (n_V(temp_Neu.V) - temp_Neu.n) / tau_nV(temp_Neu.V);

    return Neuron_Normal(temp_Neu.V, temp_Neu.n);
}

inline
double Neuron_Normal::A(double const &temp_V) const {
    return g_Na * pow(m_V(temp_V), 3) * (temp_V - V_Na);
}

inline
double Neuron_Normal::B(double const &temp_V) {
    g_Kdr = 6;

    return g_Kdr * (temp_V - V_k);
}

inline
double Neuron_Normal::C(double const &temp_V) {
    return g_L * (temp_V - V_L);
}

inline
Quartic_Equation Neuron_Normal::Solve_n(double const &temp_V) const {
    auto a = -B(temp_V);
    auto d = 1.0028 * A(temp_V);
    auto e = -(0.8820 * A(temp_V) + C(temp_V)) + Iapp;

    Quartic_Equation temp_qe;

    temp_qe.quartic_equation(a, 0, 0, d, e);

    return temp_qe;
}

inline
vector<double> Neuron_Normal::n_nullcline() const {
    vector<double> Vec_n;

    for(auto temp_V = -80.0; temp_V <= 60; temp_V += dt) {
        Vec_n.push_back(n_V(temp_V));
    }

    return Vec_n;
}

inline
vector<double> Neuron_Normal::V_nullcline() const {
    vector<Quartic_Equation> Vec_QE;
    vector<double> Vec_n;

    for (auto temp_V = -80.0; temp_V <= 60; temp_V += dt) {
        auto temp_qe = Solve_n(temp_V);
        Vec_QE.push_back(temp_qe);
    }

    for(auto &i : Vec_QE) {
        if(i.x3.imag() == 0) {
            Vec_n.push_back(i.x3.real());
        }
        else {
            Vec_n.push_back(0);
        }
    }

    return Vec_n;
}

inline
pair<vector<double>, vector<double>> Neuron_Normal::Phase_Orbit() const {
    vector<double> Vec_V;
    vector<double> Vec_n;
    Neuron_Normal temp_Neu;

    for(auto i = 0.0; i <= T; i += dt) {
        temp_Neu = HH_Model(temp_Neu);

        Vec_V.push_back(temp_Neu.V);
        Vec_n.push_back(temp_Neu.n);
    }

    return make_pair(Vec_V, Vec_n);
}

//The INa + IKdr-minimal model.
class Neuron_AD {
public:
    Neuron_AD() = default;
    Neuron_AD(double const &temp_v, double const &temp_n) : V(temp_v),  n(temp_n) { }
    virtual ~Neuron_AD() = default;

    void initialize(double const &, double const &);

    double n_V(double const &) const;
    double m_V(double const &) const;
    double tau_nV(double const &) const;
    Neuron_AD HH_Model(Neuron_AD &) const;

    double A(double const &) const;
    static double B(double const &);
    static double C(double const &);
    Quartic_Equation Solve_n(double const &) const;

    vector<double> n_nullcline() const;
    vector<double> V_nullcline() const;
    pair<vector<double>, vector<double>> Phase_Orbit() const;

protected:
    //The membrane potential.
    double V = -65;

    //The delayed rectified K+ current activation variables n.
    double n = 0.8;
    //The transient Na+ current inactivation variable h, linear relation with n.
    double h = 0.8820 - 1.0028 * n;
};

inline
void Neuron_AD::initialize(double const &temp_v, double const &temp_n) {
    V = temp_v;
    n = temp_n;
}

inline
double Neuron_AD::n_V(double const &temp_V) const {
    return 1 / (1 + exp((-43.1 - temp_V) / 7.3));
}

inline
double Neuron_AD::m_V(double const &temp_V) const {
    return 1 / (1 + exp((-30 - temp_V) / 9.5));
}

inline
double Neuron_AD::tau_nV(double const &temp_V) const {
    return 0.1 + 0.5 * 1 / (1 + exp((-27 - temp_V) / -15));
}

inline
Neuron_AD Neuron_AD::HH_Model(Neuron_AD &temp_Neu) const {
    g_Kdr = 10;

    temp_Neu.V += dt * (-g_L * (temp_Neu.V - V_L) - g_Na * pow(m_V(temp_Neu.V), 3) * temp_Neu.h * (temp_Neu.V - V_Na) - g_Kdr * pow(temp_Neu.n, 4) * (temp_Neu.V - V_k) + Iapp);
    temp_Neu.n += dt * (n_V(temp_Neu.V) - temp_Neu.n) / tau_nV(temp_Neu.V);

    return Neuron_AD(temp_Neu.V, temp_Neu.n);
}

inline
double Neuron_AD::A(double const &temp_V) const {
    return g_Na * pow(m_V(temp_V), 3) * (temp_V - V_Na);
}

inline
double Neuron_AD::B(double const &temp_V) {
    g_Kdr = 10;

    return g_Kdr * (temp_V - V_k);
}

inline
double Neuron_AD::C(double const &temp_V) {
    return g_L * (temp_V - V_L);
}

inline
Quartic_Equation Neuron_AD::Solve_n(double const &temp_V) const {
    auto a = -B(temp_V);
    auto d = 1.0028 * A(temp_V);
    auto e = -(0.8820 * A(temp_V) + C(temp_V)) + Iapp;

    Quartic_Equation temp_qe;

    temp_qe.quartic_equation(a, 0, 0, d, e);

    return temp_qe;
}

inline vector<double> Neuron_AD::n_nullcline() const {
    vector<double> Vec_n;

    for (auto temp_V = -80.0; temp_V <= 60; temp_V += dt) {
        Vec_n.push_back(n_V(temp_V));
    }

    return Vec_n;
}

inline
vector<double> Neuron_AD::V_nullcline() const {
    vector<Quartic_Equation> Vec_QE;
    vector<double> Vec_n;

    for (auto temp_V = -80.0; temp_V <= 60; temp_V += dt) {
        auto temp_qe = Solve_n(temp_V);
        Vec_QE.push_back(temp_qe);
    }

    for (auto &i : Vec_QE) {
        if (i.x3.imag() == 0) {
            Vec_n.push_back(i.x3.real());
        }
        else {
            Vec_n.push_back(0);
        }
    }

    return Vec_n;
}

inline
pair<vector<double>, vector<double>> Neuron_AD::Phase_Orbit() const {
    vector<double> Vec_V;
    vector<double> Vec_n;
    Neuron_AD temp_Neu;

    for (auto i = 0.0; i <= T; i += dt) {
        temp_Neu = HH_Model(temp_Neu);

        Vec_V.push_back(temp_Neu.V);
        Vec_n.push_back(temp_Neu.n);
    }

    return make_pair(Vec_V, Vec_n);
}

#endif _NEURON_H_
