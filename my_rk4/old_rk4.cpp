#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

using namespace std;

class ResultsWriter {
public:
    static string timestamp_now_for_filename()
    {
        std::time_t t = std::time(nullptr);
        std::tm tm;
#if defined(_MSC_VER)
        localtime_s(&tm, &t);
#else
        localtime_r(&t, &tm);
#endif
        char buf[64];
        std::strftime(buf, sizeof(buf), "%Y%m%d_%H%M%S", &tm);
        return string(buf);
    }

    static bool WriteResultsCSV(const string& out_path,
        const string& header,
        const vector<string>& param_names,
        const vector<double>& param_values,
        const vector<vector<double>>& final_results,
        int steps)
    {
        ofstream out(out_path.c_str());
        if (!out) return false;

        out << fixed << setprecision(8);
        out << header << ", ";
        for (int i = 0; i < param_names.size(); i++) {
            out << param_names[i] << "=" << param_values[i];
            if (i < param_names.size() - 1) out << ", ";
        }
        out << "\n";
        for (int row = 0; row < steps; row++) {
            for (int column = 0; column < final_results.size(); column++) {
                out << final_results[column][row];
                if (column < final_results.size() - 1) { out << ","; }
                else { out << "\n"; }
            }
        }
        out.close();

        cout << "Wrote results to: " << out_path << endl;
        return true;
    }

    static bool SaveResults(const string& title,
        const string& header,
        const vector<string>& param_names,
        const vector<double>& param_values,
        const vector<vector<double>>& final_results,
        const int steps,
        int argc, char** argv)
    {
        string out_path;
        if (argc >= 2) {
            out_path = argv[1];
        }
        else {
            string ts = timestamp_now_for_filename();
            out_path = title + "_" + ts + ".csv";
        }

        return WriteResultsCSV(out_path, header, param_names, param_values, final_results, steps);
    }

};

// initial conditions



template<typename F>
double RK4(F&& f, double t, vector<double> state, int variable, double dt)
{
    vector<double> s2 = state;
    vector<double> s3 = state;
    vector<double> s4 = state;

    double k1 = f(t, state);

    s2[variable] = state[variable] + dt * (k1 / 2);
    double k2 = f(t + dt / 2, s2);

    s3[variable] = state[variable] + dt * (k2 / 2);
    double k3 = f(t + dt / 2, s3);

    s4[variable] = state[variable] + dt * k3;
    double k4 = f(t + dt, s4);

    double y_next = state[variable] + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    return y_next;
}

double x_dot(double t, vector<double> x)
{
    double x_dot = x[0] + x[1];;
    return x_dot;
}

double y_dot(double t, vector<double> x)
{
    double y_dot = 2 * x[0] + x[1];
    return y_dot;
}

vector<double> twoD_linear_exact(double t, vector<double> x, double c1, double c2, vector<double> v1, vector<double> v2)
{
    vector<double> x_exact(2);

    x_exact[0] = c1 * exp((1 + sqrt(2)) * t) * v1[0] + c2 * exp((1 - sqrt(2)) * t) * v2[0];
    x_exact[1] = c1 * exp((1 + sqrt(2)) * t) * v1[1] + c2 * exp((1 - sqrt(2)) * t) * v2[1];

    return x_exact;
}

vector<vector<double>> numerical_jacobian(double t, const vector<double> x, double h, vector<double> rhs(double t, const vector<double> x, int n), int n, int m)
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

vector<double> twoD_linear(double t, vector<double> x, int n)
{
    vector<double> x_dot(n);

    x_dot[0] = x[0] + x[1];
    x_dot[1] = 2 * x[0] + x[1];

    return x_dot;
}


int main(int argc, char** argv) {

    vector<double> x0 = { 5.0, 10.0 };
    vector<double> v1 = { 1 / sqrt(2.0), 1.0 };
    vector<double> v2 = { -(1 / sqrt(2.0)), 1.0 };

    double c1 = (x0[1] + sqrt(2.0) * x0[0]) / 2.0;
    double c2 = (x0[1] - sqrt(2.0) * x0[0]) / 2.0;

    double t0 = 0.0;
    double dt = 0.001;
    double t_max = 9.0;
    double save_seconds = 0.25;

    int steps = round(t_max / dt);
    int save_interval_steps = ceil(save_seconds / dt);
    if (save_interval_steps < 1) save_interval_steps = 1;
    int saved_steps = ((steps + save_interval_steps - 1) / save_interval_steps) + 1;

    double t = t0;
    vector<double> x = x0;

    vector<double> x_exact(2);

    int print_interval = 250;

    cout << "Initial conditions: (x_0, y0) = (" << x[0] << ", " << x[1] << ")" << endl;

    for (int step = 0; step <= steps; step++)
    {
        t = step * dt;

        x[0] = RK4(x_dot, t, x, 0, dt);
        x[1] = RK4(y_dot, t, x, 1, dt);

        x_exact = twoD_linear_exact(t, x, c1, c2, v1, v2);
        double x_error = pow(fabs(x_exact[0] - x[0]), 2);
        double y_error = pow(fabs(x_exact[1] - x[1]), 2);


        if (step % print_interval == 0)
        {
            cout << t << ", " << x[0] << ", " << x[1] << endl;
        }
    }

    cout << "----------------------------------" << endl;

    
    // Calculate Jacobian
    double h = 1 * pow(10, -6);
    vector<vector<double>> J = numerical_jacobian(t, x, h, twoD_linear, 2, 2);
    cout << "Numerical Jacobian: " << endl;
    print_matrix(J, 2, 2);

    

    return 0;
}
