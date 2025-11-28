// Rossler_attractor.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

static struct params {
    double qT = 1.0;
    double kappa = 0.9;
    double R_ext = 500;
    double cT = 0.18;
    double i0 = 0.0;

    double V_ext = 0.0;
    double V_min = 0.0;
    double V_max = 0.0;
    double alpha_V = 0.0;

    double lam = 0.12;
    double tc = 30.0;

    double period;
};

static struct params;
void calculate_V_ext(double t, params& p);
vector<double> single_neuron_rhs(double t, vector<double> state, params& p, int n);
vector<double> RK4(double t, vector<double> y0, params& p, double dt, vector<double> rhs(double t, vector<double> y, params& p, int n), int n);

void main()
{
    int n = 4; // number of equations

    bool default_params = true;
    //bool triangular_sweep = false;


    // initial conditions
    double x0 = -1 + 0.00001;
    double T_prime0 = 0.0;
    double V0 = 0.0;
    double V_prime0 = 0.0;
    cerr << "Keep default intial conditions? (x0= " << x0 << ", T'0= " << T_prime0 << ", V0= " << V0 << ", V'0= " << V_prime0 << ")" << "|1=yes, 0=no|" << endl;
    cin >> default_params;
    if (!default_params) {
        cerr << "Enter parameters (x0, T'0, V0, V'0):" << endl;
        cin >> x0;
        cin >> T_prime0;
        cin >> V0;
        cin >> V_prime0;
    }
    vector<double> state = { x0, T_prime0, V0, V_prime0 };

    // initialise parameters

    params p;
    cerr << "Keep default params? (qT=" << p.qT << ", kappa=" << p.kappa << ", R_ext=" << p.R_ext << ", cT=" << p.cT << ", i0=" << p.i0 << ")" << " |1=yes, 0=no|" << endl;
    cin >> default_params;
    if (!default_params) {
        cerr << "Enter parameters (qT, kappa, R_ext, cT):" << endl;
        cin >> p.qT;
        cin >> p.kappa;
        cin >> p.R_ext;
        cin >> p.cT;
        cin >> p.i0;
    }

    cerr << "\nVoltage sweep parameters:" << endl;
    cerr << "V_min=";
    cin >> p.V_min;
    cerr << "V_max=";
    cin >> p.V_max;
    cerr << "period=";
    cin >> p.period;

    /*cerr << "1 or 2 way sweep?: ";
    cin >> triangular_sweep;
    if (triangular_sweep) {
        cerr << "period=";
        cin >> p.period;
    }*/

    double t = 0.0;
    double tmax, dt;
    cerr << "\ntmax = ";
    cin >> tmax;
    cerr << "dt = ";
    cin >> dt;

    p.alpha_V = (p.V_max - p.V_min) / tmax;

    int step = 0;
    int n_steps = round(tmax / dt);
    double seconds_interval;
    cerr << "print every ___ seconds?: ";
    cin >> seconds_interval;
    int print_interval = ceil(seconds_interval / dt);

    cout << "Parameters: V_min=" << p.V_min << ", V_max=" << p.V_max << ", period=" << p.period << ", qT=" << p.qT << ", kappa = " << p.kappa
        << ", R_ext = " << p.R_ext << ", cT = " << p.cT << ", i0 = " << p.i0
        << ", tc=" << p.tc << ", lam=" << p.lam << endl;
    cout << "t, " << "x, " << "T', " << "V, " << "V'" << ", V_ext" << endl;
    cerr << "|";
    while (true)
    {
        if (t > tmax)
        {
            cout << t;
            for (int i = 0; i < n; i++) {
                cout << ", " << state[i];
            }
            cerr << " Simulation Finished ==========|" << endl;
            break;
        }

        if ((step % (n_steps / 10)) == 0) {
            cerr << "=";
        }

        if (step % print_interval == 0) {
            cout << t;
            for (int i = 0; i < n; i++) {
                cout << ", " << state[i];
            }
            cout << ", " << p.V_ext;
            cout << endl;
        }

        calculate_V_ext(t, p);
        //cout << "p.V_ext=" << p.V_ext << endl;

        state = RK4(t, state, p, dt, single_neuron_rhs, n);
        t = t + dt;

        step++;
    }
}

vector<double> single_neuron_rhs(double t, vector<double> state, params& p, int n)
{
    double x = state[0];
    double T_prime = state[1];
    double V = state[2];
    double V_prime = state[3];

    double dUdx = -((x - 0.1) * (x - 0.1) + 0.1) - 180 * exp(-((x - 0.8) * (x - 0.8)) / 0.01) + 0.2 * sqrt(10.0) * p.i0 - 100.0 * pow(x, 99);

    double r = cosh(x / 0.12);
    double r_prime = sinh(x / 0.12) / 0.12;

    vector<double> x_dot(n);
    x_dot[0] = dUdx + 0.63 * V - p.qT * T_prime;
    x_dot[1] = p.cT * ((2 * V * V_prime * r - V * V * r_prime) / (r * r)) - p.kappa * T_prime;
    x_dot[2] = ((p.V_ext - (1.0 + p.R_ext / r) * V)) / p.tc;
    x_dot[3] = ((V * p.R_ext * r_prime) / (r * r) - (1 + p.R_ext / r) * V_prime) / p.tc;

    return x_dot;
}

void calculate_V_ext(double t, params& p)
{
    double periodic_t = fmod(t, p.period)/(p.period);
    if (periodic_t < 0.5) {
        p.V_ext = abs(p.V_min + 2*p.V_max * periodic_t);
    }
    else {
        p.V_ext = p.V_min + 2*(p.V_max - p.V_max * periodic_t);
    }
}

vector<double> RK4(double t, vector<double> y0, params& p, double dt, vector<double> rhs(double t, vector<double> y, params& p, int n), int n)
{
    vector<double> k1, k2, k3, k4;

    vector<double> s1 = y0;
    k1 = rhs(t, s1, p, n);

    vector<double> s2(n);
    for (int i = 0; i < n; i++) {
        s2[i] = y0[i] + dt * (k1[i] / 2);
    }
    k2 = rhs(t + dt / 2, s2, p, n);

    vector<double> s3(n);
    for (int i = 0; i < n; i++) {
        s3[i] = y0[i] + dt * (k2[i] / 2);

    }
    k3 = rhs(t + dt / 2, s3, p, n);

    vector<double> s4(n);
    for (int i = 0; i < n; i++) {
        s4[i] = y0[i] + dt * k3[i];

    }
    k4 = rhs(t + dt, s4, p, n);

    vector<double> y_next(n);
    for (int i = 0; i < n; i++) {
        y_next[i] = s1[i] + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }

    return y_next;
}
