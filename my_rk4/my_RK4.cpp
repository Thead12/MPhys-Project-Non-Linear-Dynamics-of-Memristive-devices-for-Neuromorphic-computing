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