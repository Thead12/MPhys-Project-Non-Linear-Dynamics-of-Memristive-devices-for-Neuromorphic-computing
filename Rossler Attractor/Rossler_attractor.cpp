// Rossler_attractor.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <cmath>



using namespace std;

static struct params {
    double a = 0.1;
    double b = 0.1;
    double c = 14.0;
};

static struct params;
vector<double> rossler_rhs(double t, vector<double> x, params &p, int n);
vector<double> RK4(double t, vector<double> y0, params &p, double dt, vector<double> rhs(double t, vector<double> y, params &p, int n), int n);

void main()
{
    int n = 3; // number of equations

    // initial conditions
    double x0, y0, z0;
    cerr << "Enter initial conditions (x0, y0, z0):" << endl;
    cin >> x0;
    cin >> y0;
    cin >> z0;
    vector<double> state = { x0, y0, z0 };
    
    // initialise parameters
    bool default_params;
    params p;
    cerr << "Keep default params? (a=" << p.a << ", b=" << p.b << ", c=" << p.c << ")" << "|1=yes, 0=no|" << endl;
    cin >> default_params;
    if (!default_params) {
        cerr << "Enter parameters (a, b, c):" << endl;
        cin >> p.a;
        cin >> p.b;
        cin >> p.c;
    }
    
    double t = 0.0;
    double tmax, dt;
    cerr << "tmax = ";
    cin >> tmax;
    cerr << "dt = ";
    cin >> dt;

    int step = 0;
    int n_steps = round(tmax / dt);
    double seconds_interval;
    cerr << "print every ___ seconds?: ";
    cin >> seconds_interval;
    int print_interval = ceil(seconds_interval/dt);

    cout << "t, " << "x, " << "y, " << "z" << endl;
    cerr << "|";
    while (true)
    {
        if (t > tmax)
        {
            cout << t;
            for (int i=0; i < n; i++) {
                cout << ", " << state[i];
            }
            cerr << " Simulation Finished ==========|" << endl;
            break;
        }

        if ((step % (n_steps / 10)) == 0 ) {
            cerr << "=";
        }

        if (step % print_interval == 0) {
            cout << t;
            for (int i = 0; i < n; i++) {
                cout << ", " << state[i];
            }
            cout << endl;
        }

        state = RK4(t, state, p, dt, rossler_rhs, n);
        t = t + dt;

        step++;
    }
}

vector<double> rossler_rhs(double t, vector<double> x, params &p, int n)
{
    vector<double> x_dot(n);
    x_dot[0] = -x[1] - x[2];
    x_dot[1] = x[0] + p.a * x[1];
    x_dot[2] = p.b + x[2] * (x[0] - p.c);

    return x_dot;
}

vector<double> RK4(double t, vector<double> y0, params &p, double dt, vector<double> rhs(double t, vector<double> y, params &p, int n), int n)
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
