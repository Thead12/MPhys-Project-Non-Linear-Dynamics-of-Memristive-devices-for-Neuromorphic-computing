// Rossler_attractor.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <limits>
#include <array>
#include <vector>
#include <cmath>

using namespace std;

static struct params {
    double qT = 1.0;
    double kappa = 0.9;
    double cT = 0.18;
    double i0 = 0.0;

    double lam = 0.12;
    double tc = 30.0;

    double V_ext = 0.0;
    double V_min = 0.0;
    double V_max = 0.0;
    double V_step = 0.0;

    double R_ext = 500;
    double R_min = 0.0;
    double R_max = 0.0;
    double R_step = 0.0;
};

template<size_t N>
using State = std::array<double, N>;

template<size_t N>
inline void single_neuron_rhs(double t, const State<N>& y, State<N>& dydt, const params& p);

template<size_t N, typename RHS>
inline void RK4_step(double& t, State<N>& y, const params& p, double dt, RHS&& rhs);

constexpr size_t D = 4;

int main(int argc, char* argv[])
{
    if (argc < 3) {
        std::cerr << "Usage: single_neuron.exe <params_file> <output_csv>\n";
        return 1;
    }
    std::string param_file = argv[1];
    std::string output_file = argv[2];

    std::ifstream infile(param_file);
    if (!infile.is_open()) {
        std::cerr << "Failed to open " << param_file << "\n";
        return 1;
    }

    std::ofstream outfile(output_file);
    if (!outfile.is_open()) {
        std::cerr << "Failed to open " << output_file << " for writing\n";
        return 1;
    }

    // initial conditions
    double x0 = -1 + 0.00001;
    double T_prime0 = 0.0;
    double V0 = 0.0;
    double V_prime0 = 0.0;
    State<D> state = { x0, T_prime0, V0, V_prime0 };

    params p;

    double t = 0.0;
    double tmax, dt;

    std::string line;
    while (std::getline(infile, line)) {
        // remove comments starting with #
        auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // skip empty lines
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string name;
        double value;
        if (iss >> name >> value) {
            // initial conditions
            if (name == "x0") x0 = value;
            else if (name == "T_prime0") T_prime0 = value;
            else if (name == "V0") V0 = value;
            else if (name == "V_prime0") V_prime0 = value;

            // params
            else if (name == "qT") p.qT = value;
            else if (name == "kappa") p.kappa = value;
            else if (name == "R_ext") p.R_ext = value;
            else if (name == "cT") p.cT = value;
            else if (name == "i0") p.i0 = value;
            else if (name == "lam") p.lam = value;
            else if (name == "tc") p.tc = value;

            // voltage ranges
            else if (name == "V_min") p.V_min = value;
            else if (name == "V_max") p.V_max = value;
            else if (name == "V_step") p.V_step = value;

            // Resistance ranges
            else if (name == "R_min") p.R_min = value;
            else if (name == "R_max") p.R_max = value;
            else if (name == "R_step") p.R_step = value;

            // Simulation time
            else if (name == "tmax") tmax = value;
            else if (name == "dt") dt = value;

            else std::cerr << "Unknown parameter: " << name << "\n";
        }
    }
    p.V_ext = p.V_min;
    p.R_ext = p.R_min;

    outfile << "Parameters: V_min=" << p.V_min << ", V_max=" << p.V_max << ", V_step=" << p.V_step
        << " | R_min=" << p.R_min << ", R_max=" << p.R_max << ", R_step=" << p.R_step
        << " | qT=" << p.qT << ", kappa=" << p.kappa
        << ", R_ext=" << p.R_ext << ", cT=" << p.cT << ", i0=" << p.i0
        << ", tc=" << p.tc << ", lam=" << p.lam
        << ", tmax=" << tmax << ", dt=" << dt << endl;

    outfile << "V_ext,R_ext,x_min,x_max,x_f,T'_f,V_f,V'_f" << endl; // write header

    int step = 0;
    int n_steps = round(tmax / dt);


    double x_min = 1.1;
    double x_max = -1.1;

    State<D> initial_R_min_state = state;
    while (p.V_min <= p.V_ext && p.V_ext <= p. V_max)
    {
        p.R_ext = p.R_min;

        t = 0.0;

        x_min = 1.1;
        x_max = -1.1;
        // calculating initial R_min state
        while (true)
        {
            
            if (t > tmax) // end
            {
                cerr << "-------------Incrementing V_ext=" << p.V_ext << "-------------" << endl;
                break;
            }

            // calculate
            RK4_step<D>(t, initial_R_min_state, p, dt, single_neuron_rhs<D>);
            x_min = min(x_min, initial_R_min_state[0]);
            x_max = max(x_max, initial_R_min_state[0]);
            t += dt;
        }

        state = initial_R_min_state;
        while (p.R_min <= p.R_ext && p.R_ext <= p.R_max)
        {
            cerr << "(" << p.V_ext << ", " << p.R_ext << ")" << endl;

            t = 0.0;

            x_min = 1.1;
            x_max = -1.1;

            while (true)
            {

                if (t > tmax) // end
                {
                    outfile << p.V_ext << "," << p.R_ext << "," << x_min << "," << x_max << ",";

                    // print final state
                    cerr << "inside R_ext while loop" << endl;
                    outfile << setprecision(15) << state[0];
                    for (int i = 1; i < D; i++) {
                        outfile << "," << state[i];
                    }                outfile << endl;
                    break;
                }

                // calculate
                RK4_step<D>(t, state, p, dt, single_neuron_rhs<D>);
                x_min = min(x_min, state[0]);
                x_max = max(x_max, state[0]);
                t += dt;
            }

            p.R_ext += p.R_step;
        }

        p.V_ext += p.V_step;
    }
    cerr << "\n|=============== Simulation Finished ===============|" << endl;

}
template<size_t N>
inline void single_neuron_rhs(double t, const State<N>& y, State<N>& dydt, const params& p)
{
    const double x = y[0];
    const double T_prime = y[1];
    const double V = y[2];
    const double V_prime = y[3];

    const double dUdx = -((x - 0.1) * (x - 0.1) + 0.1) - 180 * exp(-((x - 0.8) * (x - 0.8)) / 0.01) + 0.2 * sqrt(10.0) * p.i0 - 100.0 * pow(x, 99);

    const double r = cosh(x / 0.12);
    const double r2 = r * r;
    const double r_prime = sinh(x / 0.12) / 0.12;

    dydt[0] = dUdx + 0.63 * V - p.qT * T_prime;
    dydt[1] = p.cT * ((2.0 * V * V_prime * r - V * V * r_prime) / r2) - p.kappa * T_prime;
    dydt[2] = ((p.V_ext - (1.0 + p.R_ext / r) * V)) / p.tc;
    dydt[3] = ((V * p.R_ext * r_prime) / r2 - (1 + p.R_ext / r) * V_prime) / p.tc;
}

template<size_t N, typename RHS>
inline void RK4_step(double& t, State<N>& y, const params& p, double dt, RHS&& rhs)
{
    State<N> k1, k2, k3, k4, tmp;

    rhs(t, y, k1, p);
    for (int i = 0; i < N; i++) {
        tmp[i] = y[i] + dt * 0.5 * k1[i];
    }
    rhs(t + 0.5 * dt, tmp, k2, p);

    for (int i = 0; i < N; i++) {
        tmp[i] = y[i] + dt * 0.5 * k2[i];

    }
    rhs(t + 0.5 * dt, tmp, k3, p);

    for (int i = 0; i < N; i++) {
        tmp[i] = y[i] + dt * k3[i];

    }
    rhs(t + dt, tmp, k4, p);

    for (int i = 0; i < N; i++) {
        y[i] += (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}
