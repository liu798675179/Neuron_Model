#ifndef _NEURON_H_
#define _NEURON_H_

#include <random>
#include <ctime>
#include <cmath>
#include <vector>

std::default_random_engine e(time(nullptr));
std::uniform_real_distribution<double> u;

double T = 1000;
double const dt = 0.01;
double const Pi = 3.1415926;
double a = 0.02; //The time scale of U.
double b = 0.2; //The sensitivity of U.
double c = -65.0; //The after-spike reset value of V.
double d = 6.0; //The after-spike reset of U.
double const V_thr = 30.0;
double I = 70;

inline double Noise()  {
    return std::sqrt(std::abs(-4.0 * dt * std::log(u(e)) * std::cos(2 * Pi * u(e))));
}

class Neuron {
public:
    Neuron() = default;
    ~Neuron() = default;
    explicit Neuron(double const &temp_v) : V(temp_v) {}

    Neuron Izhikevich();
    void Condition();

    static std::vector<double> Simulation();
private:
    double V = -65;
    double U = 0.2 * V; //A recovery variable.
};

inline Neuron Neuron::Izhikevich() {
    V += dt * (0.04 * std::pow(V, 2) + 5 * V + 140 - U + I) + Noise();
    U += dt * a * (b * V - U);

    return Neuron(V);
}

inline void Neuron::Condition() {
    if (V >= V_thr) {
        V = c;
        U += d;
    }
}

inline std::vector<double> Neuron::Simulation(){
    Neuron temp_Neu;
    std::vector<double> Vec_Neu;

    for (double t = 0; t <= T; t += dt) {
        Vec_Neu.push_back(temp_Neu.Izhikevich().V);
        temp_Neu.Condition();
    }

    return Vec_Neu;
}


#endif _NEURON_H_
