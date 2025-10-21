#include <iostream>
#include <fstream>
#include <tuple>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
using namespace std;

#define x state[0]
#define T_prime state[1]
#define V state[2]
#define V_prime state[3]

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

extern double qT = 0.5; // "thermal charge"
extern double cT = 0.18; // "thermal capacitance"
extern double V_ext = 1.0; // "external voltage"
extern double R_ext = 500; // "external resistance"
extern double k = 0.9; // "thermal dissipation"
extern double i0 = 0.3; // "height of potential"

extern double dt = 0.001;
extern double t_max = 1;
extern double save_seconds = 0.1;

extern bool save_results = false;

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

double r(double position)
{
    return cosh(position / 0.12);
}

double r_prime(double position)
{
    return sinh(position/0.12)/0.12;
}

double pot(double position)
{
    return -(pow(position - 0.1, 2.0) + 0.1) - 180.0 * exp(-pow(position - 0.8, 2.0) / 0.01) + 0.2 * sqrt(10.0) * i0 - 100.0*(pow(position, 99.0));
}

double x_dot(double t, vector<double> state)
{
    return -pot(x) - 0.633 * V - qT * T_prime;
}

double T_prime_dot(double t, vector<double> state)
{
    return cT* (2.0 * V * V_prime * r(x) - pow(V, 2.0) * r_prime(x) / pow(r(x), 2.0)) - k * T_prime;
}

double V_dot(double t, vector<double> state)
{
    return (1.0 / 30.0) * (V_ext - (1 + R_ext) / r(x) * V);
}

double V_prime_dot(double t, vector<double> state)
{
    return 1.0 / 30.0 * (-(1 + R_ext / r(x)) * V_prime + (R_ext * r_prime(x)) / pow(r(x), 2) * V);
}

int main(int argc, char** argv) {
    int steps = round(t_max / dt);
    int save_interval_steps = ceil(save_seconds / dt);
    if (save_interval_steps < 1) save_interval_steps = 1;
    int saved_steps = ((steps + save_interval_steps - 1) / save_interval_steps) + 1;

	cout << "Total steps: " << steps << endl;
    cout << "Saved steps: " << saved_steps << endl;
    
    vector<double> t_vec(saved_steps);
    vector<double> x_vec(saved_steps);
    vector<double> T_prime_vec(saved_steps);
    vector<double> V_vec(saved_steps);
    vector<double> V_prime_vec(saved_steps);

    // initial conditions
    double t = 0.0;
    //double x_0 = -1 + 0.0001; // Causer of Nans!
    double x_0 = -0.2;
	double T_prime_0 = 0.0;
    double V_0 = 0.0;
	double V_prime_0 = 0.0;

    // state vector layout: [x, T', V, V']
    vector <double> state(4);
	state[0] = x_0;
	state[1] = T_prime_0;   
	state[2] = V_0; 
	state[3] = V_prime_0;

    cout << "Initial conditions: (x_0, T'_0, V_0, V'_0) = (" << x << ", " << T_prime << ", " << V << ", " << V_prime << ")" << endl;
    cout << "Parameters: ";
	cout << "qT=" << qT << ", cT=" << cT << ", V_ext=" << V_ext << ", R_ext=" << R_ext << ", k=" << k << ", i0=" << i0 << ", dt=" << dt << endl;
    cout << endl;
    cout << "| t  |   x   |   T'  |   V   |   V'  |" << endl;

    for (int step = 0; step <= steps; step++)
    {
        t = step * dt;

        if (step % save_interval_steps == 0)
        {
            int save_index = step / save_interval_steps;
            x_vec[save_index] = x;
            T_prime_vec[save_index] = T_prime;
            V_vec[save_index] = V;
            V_prime_vec[save_index] = V_prime;
            t_vec[save_index] = t;
            cout << t << ", " << x << ", " << T_prime << ", " << V << ", " << V_prime << endl;
        }
        x = RK4(x_dot, t, state, 0, dt);
        T_prime = RK4(T_prime_dot, t, state, 1, dt);
		V = RK4(V_dot, t, state, 2, dt);
        V_prime = RK4(V_prime_dot, t, state, 3, dt);
    }

    cout << "----------------------------------" << endl;
    cout << endl;
    
    if (save_results == false) {
        return 0;
	}
    cout << "------------------------saving results--------------------------" << endl;
	string title = "single_neuron_det_results";
	string header = "t, x, T_prime, V, V_prime";
	vector<string> param_names = { "qT", "cT", "V_ext", "R_ext", "k", "i0", "steps", "dt" };
	vector<double> param_values = { qT, cT, V_ext, R_ext, k, i0, (double)steps, dt };
    vector<vector<double>> final_results = vector<vector<double>>{t_vec, x_vec, T_prime_vec, V_vec, V_prime_vec };
    
    if (!ResultsWriter::SaveResults(title, header, param_names, param_values, final_results, saved_steps, argc, argv)) {
        cerr << "Failed to open output file" << endl;
    }
    
    return 0;
}

