// my_RK4_test.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include<array>

using namespace std;

#include "my_RK4.hpp"
#include "rk4.hpp"

int main();

void quadratic();
double quadratic_rhs(double t, double x);
double quadratic_exact(double t, double x, double c);

void linear();
double linear_rhs(double t, double x);
double linear_exact(double t, double x, double A);

void logistic();
double logistic_rhs(double t, double x);
double logistic_exact(double t, double x, double c);

void non_elementary();
double non_elementary_rhs(double t, double x);

// =========== Vector tests ==========
void twoD_linear();
vector<double> twoD_linear(double t, vector<double> x, int n);
vector<double> twoD_linear_exact(double t, vector<double> x, double c1, double c2, vector<double> v1, vector<double> v2);

void print_matrix(vector<vector<double>> A, int n, int m);


int main()
{
    cout << "\n";
    cout << "my RK4 test:\n";

    //quadratic();
    //linear();
    //logistic();
    
    //non_elementary();

    twoD_linear();
}


// ================================== Test 1: quadratic, x_dot = x^2 ==================================
void quadratic()
{
    cout << "\n--------------- x_prime = x^2 ---------------" << endl;

    // initial conditions
    double x0 = 1.0;
    double c = -(1/x0);

    double t0 = 0.0;
    double tmax = 0.5;
    double dt = 0.01;

    int n_steps = round(tmax / dt);
    int i = 0;
    int print_interval = 10;

    // delare vectors for storing results
    vector<double> t, x;
    t.reserve(n_steps), x.reserve(n_steps);
    t.push_back(t0), x.push_back(x0);

    // declare vector for x_exact
    vector<double> x_exact, x_error;
    x_exact.reserve(n_steps), x_error.reserve(n_steps);
    x_exact.push_back(x0), x_error.push_back(0.0);

 
    cout << "\nAt time t = " << t0 << ", x = " << x0 << ", dt = " << dt << endl;
    while(true)
    {
        if (t[i] > tmax)
        {
            cout << "-------------- Finished --------------" << endl;
            cout << "Reached tmax: " << t.back() << endl;
            cout << "x = " << x[i] << endl;
            cout << "x_exact = " << x_exact[i] << endl;
            break;
        }
        else if (isinf(x[i]))
        {
            cout << "\nEncountered inf at t = " << t[i] << endl;
            cout << "Previous value of x = " << x[i - 1] << endl;
            break;
        }
        else if (isnan(x[i]))
        {
            cout << "\nEncountered nan at t = " << t[i] << endl;
            cout << "Previous value of x = " << x[i - 1] << endl;
            break;
        }

        double t_next = t[i] + dt;
        double x_next = RK4(t[i], x[i], dt, quadratic_rhs);

        x_exact.push_back(quadratic_exact(t_next, x[i], c));
        x_error.push_back(pow(fabs(x_exact[i] - x[i]), 2));

        t.push_back(t_next);
        x.push_back(x_next);

        if (i % print_interval == 0)
        {
            cout << "\nStep:" << i << endl;
            cout << "RK4 result:" << endl;
            cout << "At time t = " << t[i] << ", x = " << x[i] << ", x_exact = " << x_exact[i] << endl;
            cout << "Error = " << x_error[i] << endl;
        }
        

        i++;
    }
    // calculate Mean Square Error
    double ms_error = mse(x, x_exact);
    cout << "MSE = " << ms_error << endl;
}

double quadratic_rhs(double t, double x)
{
	double x_dot = pow(x, 2);
    return x_dot;
}

double quadratic_exact(double t, double x, double c)
{
    double x_exact = -(1/(t+c));
    return x_exact;
}

// ================================== Test 2: linear, x_dot = x ==================================
void linear()
{
    cout << "\n--------------- x_prime = x ---------------" << endl;

    // initial conditions
    double x0 = 1.0;
    double A = 1.0;

    double t0 = 0.0;
    double tmax = 1.0;
    double dt = 0.01;

    int n_steps = round(tmax / dt);
    int i = 0;
    int print_interval = 20;

    // delare vectors for storing results
    vector<double> t, x;
    t.reserve(n_steps), x.reserve(n_steps);
    t.push_back(t0), x.push_back(x0);

    // declare pointer array for x_exact
    vector<double> x_exact, x_error;
    x_exact.reserve(n_steps), x_error.reserve(n_steps);
    x_exact.push_back(x0), x_error.push_back(0.0);


    cout << "\nAt time t = " << t0 << ", x = " << x0 << ", dt = " << dt << endl;
    while (true)
    {
        if (t[i] > tmax)
        {
            cout << "-------------- Finished --------------" << endl;
            cout << "Reached tmax: " << t.back() << endl;
            cout << "x = " << x[i] << endl;
            cout << "x_exact = " << x_exact[i] << endl;
            break;
        }
        else if (isinf(x[i]))
        {
            cout << "\nEncountered inf at t = " << t[i] << endl;
            cout << "Previous value of x = " << x[i - 1] << endl;
            break;
        }
        else if (isnan(x[i]))
        {
            cout << "\nEncountered nan at t = " << t[i] << endl;
            cout << "Previous value of x = " << x[i - 1] << endl;
            break;
        }

        double t_next = t[i] + dt;
        double x_next = RK4(t[i], x[i], dt, linear_rhs);

        x_exact.push_back(linear_exact(t_next, x[i], A));
        x_error.push_back(pow(fabs(x_exact[i] - x[i]), 2));

        t.push_back(t_next);
        x.push_back(x_next);

        if (i % print_interval == 0)
        {
            cout << "\nStep:" << i << endl;
            cout << "RK4 result:" << endl;
            cout << "At time t = " << t[i] << ", x = " << x[i] << ", x_exact = " << x_exact[i] << endl;
            cout << "Error = " << x_error[i] << endl;
        }


        i++;
    }
    // calculate Mean Square Error
    double ms_error = mse(x, x_exact);
    cout << "MSE = " << ms_error << endl;
}

double linear_rhs(double t, double x)
{
    double x_dot = x;
    return x_dot;
}

double linear_exact(double t, double x, double A)
{
    double x_exact = A*exp(t);
    return x_exact;
}

// ================================== Test 3: logistic, x_dot = x*(1-x) ==================================
void logistic()
{
    cout << "\n--------------- x_prime = x*(1-x) ---------------" << endl;

    // initial conditions
    double x0 = 2.0;
    double c = 1/x0 -1;

    double t0 = 0.0;
    double tmax = 2.0;
    double dt = 0.01;

    int n_steps = round(tmax / dt);
    int i = 0;
    int print_interval = 50;

    // delare vectors for storing results
    vector<double> t, x;
    t.reserve(n_steps), x.reserve(n_steps);
    t.push_back(t0), x.push_back(x0);

    // declare pointer array for x_exact
    vector<double> x_exact, x_error;
    x_exact.reserve(n_steps), x_error.reserve(n_steps);
    x_exact.push_back(x0), x_error.push_back(0.0);


    cout << "\nAt time t = " << t0 << ", x = " << x0 << ", dt = " << dt << endl;
    while (true)
    {
        if (t[i] > tmax)
        {
            cout << "-------------- Finished --------------" << endl;
            cout << "Reached tmax: " << t.back() << endl;
            cout << "x = " << x[i] << endl;
            cout << "x_exact = " << x_exact[i] << endl;
            break;
        }
        else if (isinf(x[i]))
        {
            cout << "\nEncountered inf at t = " << t[i] << endl;
            cout << "Previous value of x = " << x[i - 1] << endl;
            break;
        }
        else if (isnan(x[i]))
        {
            cout << "\nEncountered nan at t = " << t[i] << endl;
            cout << "Previous value of x = " << x[i - 1] << endl;
            break;
        }

        double t_next = t[i] + dt;
        double x_next = RK4(t[i], x[i], dt, logistic_rhs);

        x_exact.push_back(logistic_exact(t_next, x[i], c));
        x_error.push_back(pow(fabs(x_exact[i] - x[i]), 2));

        t.push_back(t_next);
        x.push_back(x_next);

        if (i % print_interval == 0)
        {
            cout << "\nStep:" << i << endl;
            cout << "RK4 result:" << endl;
            cout << "At time t = " << t[i] << ", x = " << x[i] << ", x_exact = " << x_exact[i] << endl;
            cout << "Error = " << x_error[i] << endl;
        }


        i++;
    }
    // calculate Mean Square Error
    double ms_error = mse(x, x_exact);
    cout << "MSE = " << ms_error << endl;
}

double logistic_rhs(double t, double x)
{
    double x_dot = x*(1-x);
    return x_dot;
}

double logistic_exact(double t, double x, double c)
{
    double x_exact = 1/(1 + c*exp(-t));
    return x_exact;
}

// ================================== Test 4: Non-anlytical, x_dot = ln(cos(x)-2) ==================================
void non_elementary()
{
    cout << "\n--------------- x_prime = x*(1-x) ---------------" << endl;

    // initial conditions
    double x0 = 1.0;

    double t0 = 0.0;
    double tmax = 5.0;
    double dt = 0.01;

    int n_steps = round(tmax / dt);
    int i = 0;
    int print_interval = 100;

    // delare vectors for storing results
    vector<double> t, x;
    t.reserve(n_steps), x.reserve(n_steps);
    t.push_back(t0), x.push_back(x0);

    // declare pointer array for x_exact
    vector<double> x_burkadt, x_error;
    x_burkadt.reserve(n_steps), x_error.reserve(n_steps);
    x_burkadt.push_back(x0), x_error.push_back(0.0);


    cout << "\nAt time t = " << t0 << ", x = " << x0 << ", dt = " << dt << endl;
    while (true)
    {
        if (t[i] > tmax)
        {
            cout << "-------------- Finished --------------" << endl;
            cout << "Reached tmax: " << t.back() << endl;
            cout << "x = " << x[i] << endl;
            cout << "x_burkadt = " << x_burkadt[i] << endl;
            break;
        }
        else if (isinf(x[i]))
        {
            cout << "\nEncountered inf at t = " << t[i] << endl;
            cout << "Previous value of x = " << x[i - 1] << endl;
            break;
        }
        else if (isnan(x[i]))
        {
            cout << "\nEncountered nan at t = " << t[i] << endl;
            cout << "Previous value of x = " << x[i - 1] << endl;
            break;
        }

        double t_next = t[i] + dt;
        // my RK4 
        double x_next = RK4(t[i], x[i], dt, non_elementary_rhs);

        // Burkadts library rk4
        double x_burkadt_next = rk4(t[i], x[i], dt, non_elementary_rhs);

        x_error.push_back(pow(fabs(x_burkadt_next - x_next), 2));

        t.push_back(t_next);
        x.push_back(x_next);
        x_burkadt.push_back(x_burkadt_next);

        if (i % print_interval == 0)
        {
            cout << "\nStep:" << i << endl;
            cout << "RK4 result:" << endl;
            cout << "At time t = " << t[i] << ", x = " << x[i] << endl;
            cout << "Error = " << x_error[i] << endl;
        }

        i++;
    }
    // calculate Mean Square Error
    double ms_error = mse(x, x_burkadt);
    cout << "MSE = " << ms_error << endl;
}

double non_elementary_rhs(double t, double x)
{
    double x_dot = log(2-cos(x));
    return x_dot;
}


// ================================== Vector Test 1: 2D linear, x_dot = x+y, y_dot = 2x+y ==================================
void twoD_linear()
{
    cout << "\n--------------- x_dot = x+y, y_dot = 2x+y ---------------" << endl;

    int n = 2;
    // initial conditions
    vector<double> x0 = {5.0, 10.0};
    vector<double> v1 = {1/sqrt(2.0), 1.0};
    vector<double> v2 = { -(1 / sqrt(2.0)), 1.0 };

    double c1 = (x0[1] + sqrt(2.0) * x0[0]) / 2.0;
    double c2 = (x0[1] - sqrt(2.0) * x0[0]) / 2.0;

    double t0 = 0.0;
    double tmax = 1.0;
    double dt = 0.001;

    int n_steps = round(tmax / dt);
    int i = 0;
    int print_interval = 250;

    // delare vectors for storing results
    vector<double> t;
    vector<vector<double>> x;
    t.reserve(n_steps), x.reserve(n_steps);
    t.push_back(t0), x.push_back(x0);

    // declare array for x_exact
    vector<vector<double>> x_exact;
    vector<vector<double>> x_error;
    x_exact.reserve(n_steps), x_error.reserve(n_steps);
    x_exact.push_back(x0), x_error.push_back({ 0.0, 0.0 });


    cout << "\nAt time t = " << t0 << ", x0 = " << x0[0] << ", y0 = " << x0[1] << ", dt = " << dt << endl;
    while (true)
    {
        if (t[i] > tmax)
        {
            cout << "-------------- Finished --------------" << endl;
            cout << "Reached tmax: " << t.back() << endl;
            cout << "x = " << x[i][0] << ", y = " << x[i][1] << endl;
            cout << "x_exact = " << x_exact[i][0] << ", y_exact = " << x_exact[i][1] << endl;
            break;
        }
        else if (isinf(x[i][0]) or isinf(x[i][1]))
        {
            cout << "\nEncountered inf at t = " << t[i] << endl;
            cout << "Previous value of x = " << x[i - 1][0] << ", y = " << x[i - 1][1] << endl;
            break;
        }
        else if (isnan(x[i][0]) or isnan(x[i][1]))
        {
            cout << "\nEncountered nan at t = " << t[i] << endl;
            cout << "Previous value of x = " << x[i - 1][0] << ", y = " << x[i - 1][1] << endl;
            break;
        }

        double t_next = t[i] + dt;
        vector<double> x_next = RK4(t[i], x[i], dt, twoD_linear, n);

        x_exact.push_back(twoD_linear_exact(t_next, x[i], c1, c2, v1, v2));

        double x_dif = pow(fabs(x_exact[i][0] - x[i][0]), 2);
        double y_dif = pow(fabs(x_exact[i][1] - x[i][1]), 2);
        x_error.push_back({x_dif, y_dif});

        t.push_back(t_next);
        x.push_back(x_next);

        if (i % print_interval == 0)
        {
            cout << "\nStep:" << i << endl;
            cout << "RK4 result:" << endl;
            cout << "At time t = " << t[i] << ", x = " << x[i][0] << ", y = " << x[i][1] << ", x_exact = " << x_exact[i][0] << ", y_exact = " << x_exact[i][1] << endl;
            cout << "x_error = " << x_error[i][0] << ", y_error = " << x_error[i][1] << endl;
            
            //cout << t[i] << ", " << x_error[i][0] << ", " << x_error[i][1] << endl;
        }

        i++;
    }
    // calculate Mean Square Error
    double ms_error = system_mse(x, x_exact, n);
    cout << "Total system MSE = " << ms_error << endl;
    
    // Calculate Jacobian
    double h = 1 * pow(10, -6);
    vector<vector<double>> J = numerical_jacobian(t.back(), x.back(), h, twoD_linear, n, n);
    cout << "Numerical Jacobian: " << endl;
    print_matrix(J, n, n);
}

vector<double> twoD_linear(double t, vector<double> x, int n)
{
    vector<double> x_dot(n);
    
    x_dot[0] = x[0] + x[1];
    x_dot[1] = 2 * x[0] + x[1];

    return x_dot;
}

vector<double> twoD_linear_exact(double t, vector<double> x, double c1, double c2, vector<double> v1, vector<double> v2)
{
    vector<double> x_exact(2);

    x_exact[0] = c1 * exp((1 + sqrt(2)) * t)*v1[0] + c2 * exp((1 - sqrt(2)) * t) * v2[0];
    x_exact[1] = c1 * exp((1 + sqrt(2)) * t) * v1[1] + c2 * exp((1 - sqrt(2)) * t) * v2[1];

    return x_exact;
}




void print_matrix(vector<vector<double>> A, int n, int m)
{
    cout << "\n" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << A[i][j] << ", ";
        }
        cout << endl;
    }
}

