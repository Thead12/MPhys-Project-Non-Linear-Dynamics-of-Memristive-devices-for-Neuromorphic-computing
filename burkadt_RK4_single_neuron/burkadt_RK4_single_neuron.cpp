# include <cmath>
# include <cstdlib>
# include <iomanip>
# include <iostream>

using namespace std;

# include "rk4.hpp"

int main();
void rk4vec_test();
double* rk4vec_test_f(double t, int n, double u[]);

//****************************************************************************80

int main()

{
 
    rk4vec_test();

    return 0;
}
//****************************************************************************80


void rk4vec_test()
{

    double kappa = 0.9;
    double R_ext = 500;
    double i0 = 0.3;
    double qT = 0.5;
    double cT = 0.18;

    double dt = 0.0001;
    int i;
    int n = 4;
    double t0;
    double t1;
    double tmax = 5.0;
    double* u0;
    double* u1;

    cout << "\n";
    cout << "RK4VEC_TEST\n";
    cout << "  RK4VEC takes a Runge Kutta step for a vector ODE.\n";

    t0 = 0.0;

    u0 = new double[n];
    u0[0] = -1.0+0.0001;
    u0[1] = 0.0;
    u0[2] = 195.0;
    u0[3] = 0.0;


    for (; ; )
    {
        //
        //  Print (T0,U0).
        //
        cout << "  " << setw(14) << t0
            << "  " << setw(14) << u0[0]
            << "  " << setw(14) << u0[1]
            << "  " << setw(14) << u0[2]
            << "  " << setw(14) << u0[3] << "\n";
        //
        //  Stop if we've exceeded TMAX.
        //
        if (tmax <= t0)
        {
            break;
        }
        //
        //  Otherwise, advance to time T1, and have RK4 estimate 
        //  the solution U1 there.
        //
        t1 = t0 + dt;
        u1 = rk4vec(t0, n, u0, dt, rk4vec_test_f);
        //
        //  Shift the data to prepare for another step.
        //
        t0 = t1;
        for (i = 0; i < n; i++)
        {
            u0[i] = u1[i];
        }
        delete[] u1;
    }
    return;
}
//****************************************************************************80

double* rk4vec_test_f(double t, int n, double u[])

//****************************************************************************80
//
//  Purpose:
//
//    RK4VEC_TEST_F evaluates the right hand side of a vector ODE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T, the current time.
//
//    Input, int N, the dimension of the system.
//
//    Input, double U[N], the current solution value.
//
//    Output, double RK4VEC_TEST_F[N], the value of the derivative, dU/dT.
//
{
    double kappa = 0.9;
    double R_ext = 500;
    double i0 = 0.3;
    double qT = 0.5;
    double cT = 0.18;
    double V_ext = 11.0;
    double* uprime;

    double x = u[0];
    double T_prime = u[1];
    double V = u[2];
    double V_prime = u[3];

    uprime = new double[n];

    double dUdx = -((x - 0.1) * (x - 0.1) + 0.1) - 180 * exp(-(((x - 0.8) * (x - 0.8)) / 0.01)) + 0.2 * sqrt(10) * i0 - 100 * pow(x, 99);
    double r = cosh(x / 0.12);
    double r_prime = sinh(x / 0.12) / 0.12;

    uprime[0] = dUdx-0.633*V-qT*T_prime;
    uprime[1] = cT * ((2 * V * V_prime * r - V * V * r_prime) / (r * r)) - kappa * T_prime;
    uprime[2] = (V_ext - (1 + (R_ext / r) * V)) / 30.0;
    uprime[3] = (-(1 + (R_ext / r) + (R_ext * r_prime) / (r * r))) / 30.0;

    return uprime;
}