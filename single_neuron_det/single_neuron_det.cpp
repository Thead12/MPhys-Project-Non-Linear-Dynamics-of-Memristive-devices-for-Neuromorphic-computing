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

template<typename F>
double RK4(F&& f, double t, vector<double> state, int variable, double h)
{
    vector<double> s2 = state;
    vector<double> s3 = state;
    vector<double> s4 = state;

    double k1 = f(t, state);

	s2[variable] = state[variable] + h * (k1 / 2);
    double k2 = f(t + h / 2, s2);

	s3[variable] = state[variable] + h * (k2 / 2);
    double k3 = f(t + h / 2, s3);

	s4[variable] = state[variable] + h * k3;
    double k4 = f(t + h, s4);

    double y_next = state[variable] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    return y_next;
}

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
    for (int i =0; i < param_names.size(); i++) {
        out << param_names[i] << "=" << param_values[i];
        if (i < param_names.size() - 1) out << ", ";
	}
    out << "\n";
    for (int row = 0; row < steps; row++) {
        for (int column = 0; column < final_results.size(); column++) {
            out << final_results[column][row];
            if (column < final_results.size() - 1) {out << ",";} else {out << "\n";}
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

extern double qT = 0.5; // "thermal charge"
extern double cT = 0.18; // "thermal capacitance"
extern double V_ext = 0.0; // "external voltage"
extern double R_ext = 500; // "external resistance"
extern double k = 0.9; // "thermal dissipation"
extern double i0 = 0.3; // "height of potential"

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
    int steps = 1000;
    double h = 0.01;
    double t = 0.0;

    vector<double> t_vec(steps);
    vector<double> x_vec(steps);
    vector<double> T_prime_vec(steps);
    vector<double> V_vec(steps);
    vector<double> V_prime_vec(steps);

    // initial conditions
    double x_0 = -1 + 0.0001;
	double T_prime_0 = 0.0;
    double V_0 = 0.0;
	double V_prime_0 = 0.0;

    // state vector layout: [x, T', V, V']
    vector <double> state(4);
	state[0] = x_0;
	state[1] = T_prime_0;   
	state[2] = V_0; 
	state[3] = V_prime_0;

    cout << "Initial conditions: (x, T', V, V') = (" << x << ", " << T_prime << ", " << V << ", " << V_prime << ")" << endl;
    cout << "Parameters: ";
	cout << "qT=" << qT << ", cT=" << cT << ", V_ext=" << V_ext << ", R_ext=" << R_ext << ", k=" << k << ", i0=" << i0 << endl;
    cout << endl;
    cout << "| t  |   x   |   T'  |   V   |   V'  |" << endl;
    for (int n = 0; n < steps; n++)
    {
        x = RK4(x_dot, t, state, 0, h);
        T_prime = RK4(T_prime_dot, t, state, 1, h);
		V = RK4(V_dot, t, state, 2, h);
        V_prime = RK4(V_prime_dot, t, state, 3, h);

        x_vec[n] = x;
		T_prime_vec[n] = T_prime;
		V_vec[n] = V;
		V_prime_vec[n] = V_prime;

        V_ext += 1*h;

        if (n % 10 == 0) {
            cout << t << ", " << x << ", " << T_prime << ", " << V << ", " << V_prime << ", V_ext = " << V_ext << endl;
        }
		t += h;
		t_vec[n] = t;
    }

    cout << "----------------------------------" << endl;
    cout << endl;
    
	// ------------------------saving results--------------------------
	string title = "single_neuron_det_results";
	string header = "t, x, T_prime, V, V_prime";
	vector<string> param_names = { "qT", "cT", "V_ext", "R_ext", "k", "i0", "steps", "h" };
	vector<double> param_values = { qT, cT, V_ext, R_ext, k, i0, (double)steps, h };
    vector<vector<double>> final_results = vector<vector<double>>{t_vec, x_vec, T_prime_vec, V_vec, V_prime_vec };
    
    if (!SaveResults(title, header, param_names, param_values, final_results, steps, argc, argv)) {
        cerr << "Failed to open output file" << endl;
    }
    
    return 0;
}
