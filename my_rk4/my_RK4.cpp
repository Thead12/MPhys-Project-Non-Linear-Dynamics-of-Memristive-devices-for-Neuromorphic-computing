#include <cmath>
#include <vector>

using namespace std;

#include "my_RK4.hpp"

double RK4(double t, double y0, double dt, double rhs(double t, double y))
{
    double s1, s2, s3, s4;
    double k1, k2, k3, k4;

    s1 = y0;
    k1 = rhs(t, s1);

    s2 = y0 + dt * (k1 / 2);
    k2 = rhs(t + dt / 2, s2);

    s3 = y0 + dt * (k2 / 2);
    k3 = rhs(t + dt / 2, s3);

    s4 = y0 + dt * k3;
    k4 = rhs(t + dt, s4);

    double y_next = s1 + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

    return y_next;
}

vector<double> RK4(double t, vector<double> y0, double dt, vector<double> rhs(double t, vector<double> y, int n), int n)
{
    vector<double> k1, k2, k3, k4;

    vector<double> s1 = y0;
    k1 = rhs(t, s1, n);

    vector<double> s2(n);
    for (int i = 0; i < n; i++) {
        s2[i] = y0[i] + dt * (k1[i] / 2);
    }
    k2 = rhs(t + dt / 2, s2, n);

    vector<double> s3(n);
    for (int i = 0; i < n; i++) {
        s3[i] = y0[i] + dt * (k2[i] / 2);

    }
    k3 = rhs(t + dt / 2, s3, n);

    vector<double> s4(n);
    for (int i = 0; i < n; i++) {
        s4[i] = y0[i] + dt * k3[i];

    }
    k4 = rhs(t + dt, s4, n);

    vector<double> y_next(n);
    for (int i = 0; i < n; i++) {
        y_next[i] = s1[i] + (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }

    return y_next;
}

vector<vector<double>> numerical_jacobian(double t, const vector<double> x, double h,vector<double> rhs(double t, const vector<double> x, int n), int n, int m)
{
    vector<vector<double>> J(m, vector<double>(n));
    vector<double> f = rhs(t, x, n);

    for (int j = 0; j < n; j++)  
    {
        vector<double> x_pert = x;
        x_pert[j] += h;

        vector<double> f_pert = rhs(t, x_pert, n);

        for (int i = 0; i < m; i++)  
        {
            J[i][j] = (f_pert[i] - f[i]) / h; 
        }
    }

    return J;
}



double mse(vector<double> x, vector<double> y)
{
    if (x.empty() or y.empty())
    {
        return -1.0;
    }
    else if (x.size() != y.size())
    {
        return -2.0;
    }

    double sum = 0.0;
    for (int i = 0; i < x.size(); i++)
    {
        double error = fabs(y[i] - x[i]);
        sum = pow(error, 2);
    }
    double mean = sum / x.size();
    return mean;
}

double system_mse(vector<vector<double>> x, vector<vector<double>> y, int n)
{
    if (x.empty() or y.empty())
    {
        return -1.0;
    }
    else if (x.size() != y.size())
    {
        return -2.0;
    }

    vector<double> x_mse(n);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < x.size(); j++)
        {
            double error = fabs(y[j][i] - x[j][i]);
            double temp = pow(error, 2);
            x_mse[i] += temp;
        }
        x_mse[i] = x_mse[i] / x.size();   
    }

    double sum = 0.0;
    for (int i = 0; i < n; i++)
    {
        sum += x_mse[i];
    }

    double total_mse = sum / n;
    return total_mse;
}




