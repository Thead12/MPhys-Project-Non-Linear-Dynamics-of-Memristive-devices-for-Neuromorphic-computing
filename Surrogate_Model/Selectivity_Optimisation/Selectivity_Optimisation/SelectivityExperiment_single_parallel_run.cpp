#include <iostream>
#include <string>
#include <limits>
#include <array>
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>

//# include <pybind11/pybind11.h>
//# include <pybind11/numpy.h>

# define M_PI           3.14159265358979323846  /* pi */

//namespace py = pybind11;

using namespace std;

// Set DEBUG_OUTPUT to true to enable sweep summary and selectivity debug prints.
// The progress bar always remains enabled.
#ifdef NDEBUG
constexpr bool DEBUG_OUTPUT = false;
#else
constexpr bool DEBUG_OUTPUT = true;
#endif

struct params {
    double qT = 1.0;
    double kappa = 0.9;
    double cT = 0.18;
    double i0 = 0.0;

    double lam = 0.12;
    double tc = 30.0;

    double V_ext = 0.0;
    double V_AC;

    double V_DC = 6.0;
    double R_ext = 500.0;

    double warmUp_cycles; 
    double num_cycles;     
    double dt;

    double target_freq;
    double frequency_range;
    double frequency_step;

    double depth_checking_range;

    double spike_threshold = -0.8;

};

template<size_t N>
using State = std::array<double, N>;

template<size_t N>
inline void single_neuron_rhs(double t, const State<N>& y, State<N>& dydt, const params& p, double V_ext_val);

template<size_t N, typename RHS>
inline void RK4_step(double& t, State<N>& y, const params& p, double dt, RHS&& rhs, double V_ext_val);

constexpr size_t D = 4;

// Main simulation function that returns observable measurements
struct SimulationObservables {
    double width;
    double depth;
    double CV1_at_target_freq;
    double CV2_at_target_freq;
    double centre_of_phase_locking;
};

struct ParametricSweepResult {
    double firing_rate;
    double isi_mean;
    double winding_number;
    double driving_frequency;
    double CV1;
    double CV2;
};

double calculateSD(const std::vector<double>& data, double& mean) {
    double sum = 0.0;
    double standardDeviation = 0.0;

    for (double value : data) {
        sum += value;
    }
    mean = sum / data.size();

    for (double value : data) {
        standardDeviation += pow(value - mean, 2);
    }

    return sqrt(standardDeviation / (data.size() - 1));
}

double calculateCV2(const std::vector<double>& t) {
    double sum_val = 0.0;
    for (int n = 0; n < t.size() - 1; n++)
    {
        double abs_err = abs(t[n + 1] - t[n]);
        double tot_interval = t[n + 1] + t[n];
        sum_val += 2.0 * (abs_err / tot_interval);
    }
    return sum_val / (t.size() - 1);

}

//py::dict run_selectivity_simulation(double V_DC, double R_ext, double V_AC = 1.0,
void run_selectivity_simulation(double V_DC, double R_ext, double V_AC=3.0,
    double target_freq=0.064,
    double frequency_range=0.04,
    double frequency_step=0.001,
    double warmUp_cycles=10.0,
    double num_cycles=100.0,
    double depth_checking_range=0.03,
    double dt=1e-3)
{
    cerr << target_freq << endl;
    // initial conditions
    double x0 = -1 + 0.00001;
    double T_prime0 = 0.0;
    double V0 = 0.0;
    double V_prime0 = 0.0;
    State<D> initial_state = { x0, T_prime0, V0, V_prime0 };

    params p;
    p.V_DC = V_DC;
    p.R_ext = R_ext;
    p.V_AC = V_AC;
    p.target_freq = target_freq;
    p.frequency_range = frequency_range;
    p.frequency_step = frequency_step;
    p.warmUp_cycles = warmUp_cycles;
    p.num_cycles = num_cycles;
    p.depth_checking_range = depth_checking_range;
    p.dt = dt;

    

    double f_min = std::clamp(p.target_freq - p.frequency_range, 0.001, 1e5);
    double f_max = std::clamp(p.target_freq + p.frequency_range, 0.001, 1e5);

    int f_steps = static_cast<int>(std::round(std::abs(f_max - f_min) / std::abs(p.frequency_step)));

    std::vector<ParametricSweepResult> sweep_results(f_steps + 1);

    // frequency iterations
#pragma omp parallel for
    for (int i = 0; i <= f_steps; ++i)
    {
        if constexpr (DEBUG_OUTPUT)
        {
            if (i == 0) { // Only print once to avoid flooding the console
                int total_threads = omp_get_num_threads();
                std::cout << "Running on " << total_threads << " threads!" << std::endl;
            }
        }
        State<D> current_state = { x0, T_prime0, V0, V_prime0 };

        double x_min = 1.1;
        double x_max = -1.1;

        bool have_we_spiked = false;
        int spike_count = 0;
        double t = 0.0;
        double t_warmUp = 0.0;

        std::vector<double> spike_timings;
        spike_timings.reserve(static_cast<int>(p.num_cycles * 4));

        double driving_frequency = f_min + (i * p.frequency_step);
        double angular_frequency = 2 * M_PI * driving_frequency;

        double tmax = (p.num_cycles) / (driving_frequency);
        int t_steps = static_cast<int>(std::round(tmax / p.dt));
        int warmUp_steps = static_cast<int>(std::round((p.warmUp_cycles / driving_frequency) / p.dt));

        t = 0.0;
        x_min = 100.0;
        x_max = -100.0;
        spike_count = 0;
        have_we_spiked = false;

        // reset to initial state variables
        current_state = initial_state;
        // warm up without taking spiking measurements
        double V_ext_val = 0.0;
        for (int warmUp_step = 0; warmUp_step < warmUp_steps; ++warmUp_step)
        {
            t_warmUp = warmUp_step * p.dt;
            V_ext_val = p.V_DC + p.V_AC * cos(angular_frequency * t);
            RK4_step<D>(t, current_state, p, p.dt, single_neuron_rhs<D>, V_ext_val);
        }

        // Begin time stepping - collect data
        for (int k = 0; k < t_steps; ++k)
        {
            t = t_warmUp + (k * p.dt);
            V_ext_val = p.V_DC + p.V_AC * cos(angular_frequency * t);
            RK4_step<D>(t, current_state, p, p.dt, single_neuron_rhs<D>, V_ext_val);

            if (current_state[0] > p.spike_threshold && !have_we_spiked) {
                spike_count++;
                spike_timings.push_back(t);
                have_we_spiked = true;
            }
            else if (current_state[0] < p.spike_threshold && have_we_spiked) {
                have_we_spiked = false;
            }

            x_min = min(x_min, current_state[0]);
            x_max = max(x_max, current_state[0]);
        }

        std::vector<double> isi(size(spike_timings) - 1);

        for (int i = 0; i < size(spike_timings) - 1; ++i)
        {
            isi[i] = spike_timings[i + 1] - spike_timings[i];
        }

        ParametricSweepResult sweep_result;

        double isi_mean = 0.0;
        double isi_standard_deviation = 0.0;
        double CV1 = 0.0;
        double CV2 = 0.0;

        if (spike_timings.size() >= 2) {
            std::vector<double> isi(spike_timings.size() - 1);
            for (size_t i = 0; i < isi.size(); ++i) {
                isi[i] = spike_timings[i + 1] - spike_timings[i];
            }
            isi_standard_deviation = calculateSD(isi, isi_mean);
            if (isi_mean > 0.0) {
                CV1 = isi_standard_deviation / isi_mean;
                CV2 = calculateCV2(isi);
            }
        }

        sweep_result.winding_number = spike_count / p.num_cycles;
        sweep_result.driving_frequency = driving_frequency;
        sweep_result.CV1 = CV1;
        sweep_result.CV2 = CV2;
        sweep_result.firing_rate = spike_count / tmax;
        sweep_result.isi_mean = isi_mean;

        sweep_results[i] = sweep_result;

        t = t_steps * p.dt;
    }

    if constexpr (DEBUG_OUTPUT) {
        cerr << "\n|=============== Simulation Finished ===============|" << endl;
        cerr << "V_DC= " << p.V_DC << ", R_ext= " << p.R_ext << ", Driving Frequency= " << p.target_freq << endl;
        cerr << "driving_freq,firing_rate,ISI_mean,CV1,CV2,winding_number" << endl;
        for (const auto& result : sweep_results)
        {
            cerr << result.driving_frequency << "," << result.firing_rate << "," << result.isi_mean << "," << result.CV1 << "," << result.CV2 << "," << result.winding_number << endl;
        }
    }

    // calculating width and depth of 1:1 phase locking region
    double width = 0.0;
    std::vector<double> locking_frequencies;

    locking_frequencies.reserve(f_steps);

    std::vector<int> locking_indices;
    locking_indices.reserve(f_steps);
    for (int i = 0; i < f_steps; i++)
    {
        const auto& result = sweep_results[i];
        if (abs(result.winding_number - 1.0) < 0.05) // Check for 1:1 phase locking
        {
            width += p.frequency_step; // Increment width by the frequency step size
            locking_frequencies.push_back(result.driving_frequency);
            locking_indices.push_back(i);
        }
    }

    int is_phase_lock = 1;
    double start_locking_freq = 0.0;
    double end_locking_freq = 0.0;
    double phase_locking_centre = 0.0;
    double left_depth = 0.0;
    double right_depth = 0.0;

    if (!locking_frequencies.empty()) {
        start_locking_freq = locking_frequencies.front();
        end_locking_freq = locking_frequencies.back();
        phase_locking_centre = (start_locking_freq + end_locking_freq) / 2.0;

        int max_index = static_cast<int>(sweep_results.size());

        int left_edge_index = locking_indices.front();
        int right_edge_index = locking_indices.back();

        double left_edge_freq = locking_frequencies.front();
        double right_edge_freq = locking_frequencies.back();

        // Calculate target frequency for left depth measurement
        double left_target_freq = left_edge_freq - p.depth_checking_range;
        // Find index offset based on frequency difference

        int left_index_offset = static_cast<int>(std::round((left_target_freq - left_edge_freq) / p.frequency_step));
        int left_depth_index = std::clamp(left_edge_index + left_index_offset, 0, max_index);

        // Bounds checking for left depth
        if (left_depth_index >= 0 && left_depth_index < max_index) {
            left_depth = sweep_results[left_depth_index].CV2;
        }
        else if (left_depth_index < 0) {
            left_depth = sweep_results[0].CV2; // Use first available
        }
        else {
            left_depth = sweep_results.back().CV2; // Use last available
        }

        // Calculate target frequency for right depth measurement
        double right_target_freq = right_edge_freq + p.depth_checking_range;
        // Find index offset based on frequency difference
        int right_index_offset = static_cast<int>(std::round((right_target_freq - right_edge_freq) / p.frequency_step));
        int right_depth_index = std::clamp(right_edge_index + right_index_offset, 0, max_index);

        // Bounds checking for right depth
        if (right_depth_index >= 0 && right_depth_index < max_index) {
            right_depth = sweep_results[right_depth_index].CV2;
        }
        else if (right_depth_index < 0) {
            right_depth = sweep_results[0].CV2; // Use first available
        }
        else {
            right_depth = sweep_results.back().CV2; // Use last available
        }

        if constexpr (DEBUG_OUTPUT) {
            cout << "1:1 phase locking region width: " << width << " Hz" << endl;
            cout << "1:1 phase locking region start frequency: " << start_locking_freq << " Hz" << endl;
            cout << "1:1 phase locking region end frequency: " << end_locking_freq << " Hz" << endl;
            cout << "1:1 phase locking region centre frequency: " << phase_locking_centre << " Hz" << endl;
            cout << "index of left depth measurement: " << left_depth_index << " (frequency: " << sweep_results[left_depth_index].driving_frequency << " Hz)" << endl;
            cout << "index of right depth measurement: " << right_depth_index << " (frequency: " << sweep_results[right_depth_index].driving_frequency << " Hz)" << endl;
            cout << "Depth on left side of locking region: " << left_depth << endl;
            cout << "Depth on right side of locking region: " << right_depth << endl;
            cout << "Average depth (for PPO reward): " << (left_depth + right_depth) / 2.0 << endl;
        }

    }
    else {
        width = 999.0;
        left_depth = 999.0;
        right_depth = 999.0;
        phase_locking_centre = 999.0;
        is_phase_lock = 0;
        if constexpr (DEBUG_OUTPUT) {
            cout << "No 1:1 phase locking region found" << endl;
        }
    }

    double CV1_at_target_freq = sweep_results[f_steps / 2].CV1;
    double CV2_at_target_freq = sweep_results[f_steps / 2].CV2;

    /*
    py::dict result;
    result["is_phase_lock"] = is_phase_lock;
    result["width"] = width;
    result["CV1_at_target_freq"] = CV1_at_target_freq;
    result["CV2_at_target_freq"] = CV2_at_target_freq;
    result["centre_of_phase_locking"] = phase_locking_centre;
    result["left_depth"] = left_depth;
    result["right_depth"] = right_depth;
    result["depth"] = (left_depth + right_depth) / 2.0;  // Average depth for PPO reward function
    return result;
    */
}

template<size_t N>
inline void single_neuron_rhs(double t, const State<N>& y, State<N>& dydt, const params& p, double V_ext_val)
{
    const double x = y[0];
    const double T_prime = y[1];
    const double V = y[2];
    const double V_prime = y[3];

    double x2 = x * x;
    double x4 = x2 * x2;
    double x8 = x4 * x4;
    double x16 = x8 * x8;

    const double dUdx = -((x - 0.1) * (x - 0.1) + 0.1) - 180 * exp(-((x - 0.8) * (x - 0.8)) / 0.01) + 0.2 * sqrt(10.0) * p.i0 - 100 * pow(x, 99);

    const double inv_012 = 1.0 / 0.12;
    const double arg = x * inv_012;
    const double exp_pos = std::exp(arg);
    const double exp_neg = 1.0 / exp_pos;

    const double r = 0.5 * (exp_pos + exp_neg); // replaces cosh(x / 0.12)
    const double r_inv = 1.0 / r;
    const double r_inv2 = r_inv * r_inv;
    const double r_prime = 0.5 * (exp_pos - exp_neg) * inv_012; // replaces sinh(x / 0.12) / 0.12

    dydt[0] = dUdx + 0.63 * V - p.qT * T_prime;
    dydt[1] = p.cT * ((2.0 * V * V_prime * r - V * V * r_prime) * r_inv2) - p.kappa * T_prime;
    dydt[2] = ((V_ext_val - (1.0 + (p.R_ext * r_inv)) * V)) / p.tc;
    dydt[3] = (((V * p.R_ext * r_prime) * r_inv2) - (1 + (p.R_ext * r_inv)) * V_prime) / p.tc;
}

template<size_t N, typename RHS>
inline void RK4_step(double& t, State<N>& y, const params& p, double dt, RHS&& rhs, double V_ext_val)
{
    State<N> k1, k2, k3, k4, tmp;

    rhs(t, y, k1, p, V_ext_val);
    tmp[0] = y[0] + dt * 0.5 * k1[0];
    tmp[1] = y[1] + dt * 0.5 * k1[1];
    tmp[2] = y[2] + dt * 0.5 * k1[2];
    tmp[3] = y[3] + dt * 0.5 * k1[3];

    rhs(t + 0.5 * dt, tmp, k2, p, V_ext_val);
    tmp[0] = y[0] + dt * 0.5 * k2[0];
    tmp[1] = y[1] + dt * 0.5 * k2[1];
    tmp[2] = y[2] + dt * 0.5 * k2[2];
    tmp[3] = y[3] + dt * 0.5 * k2[3];

    rhs(t + 0.5 * dt, tmp, k3, p, V_ext_val);
    tmp[0] = y[0] + dt * k3[0];
    tmp[1] = y[1] + dt * k3[1];
    tmp[2] = y[2] + dt * k3[2];
    tmp[3] = y[3] + dt * k3[3];

    rhs(t + dt, tmp, k4, p, V_ext_val);
    y[0] += (dt / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    y[1] += (dt / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    y[2] += (dt / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    y[3] += (dt / 6.0) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
}

/*
template<size_t N, typename RHS>
inline void RK4_step(double& t, State<N>& y, const params& p, double dt, RHS&& rhs)
{
    State<N> k1, k2, k3, k4, tmp;

    rhs(t, y, k1, p);
    for (int i = 0; i < N; i++) {
        tmp[i] = y[i] + dt * 0.5 * k1[i];
    }
    rhs(t + 0.5 * dt, tmp, k2, p);

    for (int i = 0; i < N; i++) {
        tmp[i] = y[i] + dt * 0.5 * k2[i];
    }
    rhs(t + 0.5 * dt, tmp, k3, p);

    for (int i = 0; i < N; i++) {
        tmp[i] = y[i] + dt * k3[i];
    }
    rhs(t + dt, tmp, k4, p);

    for (int i = 0; i < N; i++) {
        y[i] += (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    }
}
*/

int main() {
    std::cout << "Starting Debug Run..." << std::endl;

    // Call your function with sample parameters
    // Note: I'm using sample values based on your previous code
    run_selectivity_simulation(6.0, 500.0, 1.0);

    // Print something so you can see if it worked
    std::cout << "Simulation Finished!" << std::endl;

    return 0;
}

/*
PYBIND11_MODULE(SelectivityExperiment, m) {
    m.doc() = "Neuron selectivity simulation using RK4 integration";
    m.def("run_simulation", &run_selectivity_simulation,
        py::arg("V_DC"), py::arg("R_ext"), py::arg("V_AC") = 1.0,
        py::arg("target_freq") = 0.064, py::arg("frequency_range") = 0.02, py::arg("frequency_step") = 0.001,
        py::arg("warmUp_cycles") = 10.0, py::arg("num_cycles") = 50.0, py::arg("depth_checking_range") = 0.002, py::arg("dt") = 1e-3,
        "Run selectivity simulation and return macroscopic measurements.\n"
        "Args:\n"
        "  V_DC: DC voltage\n"
        "  R_ext: External resistance\n"
        "  V_AC: AC voltage\n"
        "  target_freq: Target driving frequency (Hz)\n"
        "  frequency_range: Range around target frequency to sweep (Hz)\n"
        "  frequency_step: Step size for frequency sweep (Hz)\n"
        "  warmUp_cycles: Number of cycles to run for warm-up (no measurements)\n"
        "  num_cycles: Number of cycles to run for measurements\n"
        "  depth_checking_range: Range for checking depth\n"
        "  dt: Time step for integration\n"
        "Returns:\n"
        "  dict with keys: width, CV1_at_target_freq, CV2_at_target_freq, centre_of_phase_locking, left_depth, right_depth, depth");
}
*/