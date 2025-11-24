#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

using namespace std;

#include "my_RK4.hpp"

int single_neuron_det(int argc, char** argv);
vector<double> single_neuron_det_rhs(double t, vector<double> x, int n);

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

        for (int row = 0; row < final_results.size(); row++) {
            for (int column = 0; column < final_results[0].size(); column++) {
                out << final_results[row][column];
                if (column < final_results[0].size() - 1) { out << ","; }
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

void print_matrix(vector<vector<double>> A, int n, int m)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << A[i][j] << ", ";
        }
        cout << endl;
    }
}

int main(int argc, char** argv) {
    single_neuron_det(argc, argv);
}

// constants
const double qT = 0.5; // "thermal charge"
const double cT = 0.18; // "thermal capacitance"
const double V_ext = 1.8; // "external voltage"
const double R_ext = 500; // "external resistance"
const double kappa = 0.9; // "thermal dissipation"
const double i0 = 0.3; // "height of potential"

int single_neuron_det(int argc, char** argv)
{
    int n = 4;

    // initial conditions
    //double x_0 = -1 + 0.00001;
    double x_0 = -0.85;

    double T_prime_0 = 0.0;
    double V_0 = 30.0;
    double V_prime_0 = 0.0;

    double t0 = 0.0;
    double t_max = 2.0;
    double dt = 1.0*pow(10, -6);
    double save_seconds = 0.01;

    int n_steps = round(t_max / dt);
    int save_interval_steps = ceil(save_seconds / dt);
    if (save_interval_steps < 1) save_interval_steps = 1;
    int saved_steps = ((n_steps + save_interval_steps - 1) / save_interval_steps) + 1;
    bool save_results = true;

    vector<double> t;
    vector<vector<double>> x;
    t.reserve(saved_steps);
    x.reserve(saved_steps);
    t.push_back(t0), x.push_back({ x_0, T_prime_0, V_0, V_prime_0 });

    double t_next = 0.0;
    vector<double> x_next;
    
    double t_curr = 0.0;
    vector<double> x_curr = x.back();

    cout << "Initial conditions: (x_0, T'_0, V_0, V'_0) = (" << x_0 << ", " << T_prime_0 << ", " << V_0 << ", " << V_prime_0 << ")" << endl;
    cout << "Parameters: ";
    cout << "qT=" << qT << ", cT=" << cT << ", V_ext=" << V_ext << ", R_ext=" << R_ext << ", kappa=" << kappa << ", i0=" << i0 << ", dt=" << dt << endl;
    cout << endl;
    cout << "| t  |   x   |   T'  |   V   |   V'  |" << endl;
   
    int i = 0;
    while (true)
    {
        if (t_next > t_max)
        {
            cout << "-------------- Finished --------------" << endl;
            cout << "Reached tmax: " << t_next << endl;
            cout << t_next << ", " << x.back()[0] << ", " << x.back()[1] << ", " << x.back()[2] << ", " << x.back()[3] << endl;

            break;
        }


        x_next = RK4(t_curr, x_curr, dt, single_neuron_det_rhs, n);
        t_next = t_curr + dt;

        x_curr = x_next;
        t_curr = t_next;

        if (i % save_interval_steps == 0)
        {
            cout << t.back() << ", " << x.back()[0] << ", " << x.back()[1] << ", " << x.back()[2] << ", " << x.back()[3] << endl;

            t.push_back(t_next);
            x.push_back(x_next);

        }
        i++;
    }

    // Calculate Jacobian
    double h = 1 * pow(10, -6);
    vector<vector<double>> J = numerical_jacobian(t.back(), x.back(), h, single_neuron_det_rhs, n, n);
    cout << "\nNumerical Jacobian: " << endl;
    print_matrix(J, n, n);

    if (save_results == false) {
        return 0;
    }
    cout << "------------------------saving results--------------------------" << endl;
    string title = "single_neuron_det_results";
    string header = "t, x, T_prime, V, V_prime";
    vector<string> param_names = { "qT", "cT", "V_ext", "R_ext", "k", "i0", "dt", "x0", "T'0", "V0", "V'0" };
    vector<double> param_values = { qT, cT, V_ext, R_ext, kappa, i0, dt, x_0, T_prime_0, V_0, V_prime_0 };

    vector<vector<double>> final_results = x;

    for (int i = 0; i < final_results.size(); i++)
    {
        final_results[i].insert(final_results[i].begin(), t[i]);
    }
    
    cout << final_results.size() << ", " << final_results[0].size() << ", " << final_results[1].size() << ", " << final_results.back().size() << endl;;

    if (!ResultsWriter::SaveResults(title, header, param_names, param_values, final_results, saved_steps, argc, argv)) {
        cerr << "Failed to open output file" << endl;
    }

    return 0;

    
}

vector<double> single_neuron_det_rhs(double t, vector<double> state, int n)
{
    double x = state[0];
    double T_prime = state[1];
    double V = state[2];
    double V_prime = state[3];

    double dUdx = -(pow(x - 0.1, 2.0) + 0.1) - 180.0 * exp(-pow(x - 0.8, 2.0) / 0.01) + 0.2 * sqrt(10.0) * i0 - 100.0 * (pow(x, 99.0));

    double r = cosh(x / 0.12);
    double r_prime = sinh(x / 0.12) / 0.12;

    vector<double> x_dot(n);
    x_dot[0] = dUdx - 0.633 * V - qT * T_prime;
    x_dot[1] = cT * ((2.0 * V * V_prime * r - (V*V) * r_prime) / (r*r)) - kappa * T_prime;
    x_dot[2] = (V_ext - (1.0 + (R_ext / r)) * V)/30.0;
    x_dot[3] = (-(1.0 + (R_ext / r)) * V_prime + ((R_ext * r_prime) / (r*r)) * V)/30.0;

    return x_dot;
}


