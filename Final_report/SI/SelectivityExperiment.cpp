#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <omp.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h> 

#define M_PI 3.14159265358979323846

using namespace std;
namespace py = pybind11;

#ifdef NDEBUG
constexpr bool DEBUG_OUTPUT = false;
#else
constexpr bool DEBUG_OUTPUT = true;
#endif

/*
C++ numerical continuation, integration, and selectivity analysis code
*/

constexpr size_t D = 4;

struct params {
    double qT = 1.0;
    double kappa = 0.9;
    double cT = 0.18;
    double i0 = 0.0;

    double lam = 0.12;
    double tc = 30.0;

    double V_ext = 0.0;
    double V_AC = 1.0;

    double V_DC = 6.0;
    double R_ext = 500.0;

    double warmUp_spikes = 100.0;
    double max_spikes = 200.0;
    double dt = 1e-3;

    double target_freq = 0.06;
    double frequency_range = 0.04;
    double frequency_step = 0.0005;

    double depth_checking_range = 0.0025;

    double spike_threshold = -0.8;
};

struct SimulationInput {
    double V_DC=-1.0;
    double R_ext=-1.0;
};

struct SimulationOutput {
    int is_phase_lock=-1;
    double frequency=-1.0;
    double width=-1.0;
    double CV1_at_target_freq=-1.0;
    double CV2_at_target_freq=-1.0;
    double centre_of_phase_locking=-1.0;
    double left_depth=-1.0;
    double right_depth=-1.0;
    double depth=-1.0;
};

struct ParametricSweepResult {
    int spike_count=-1.0;
    double tmax = -1.0;
    double firing_rate=-1.0;
    double isi_mean=-1.0;
    double winding_number=-1.0;
    double driving_frequency=-1.0;
    double CV1=-1.0;
    double CV2=-1.0;
};

template<size_t N>
using State = std::array<double, N>;

double calculateSD(const std::vector<double>& data, double& mean) {
    double sum = 0.0;
    double standardDeviation = 0.0;
    for (double value : data) sum += value;
    mean = sum / data.size();
    for (double value : data) standardDeviation += pow(value - mean, 2);
    return sqrt(standardDeviation / (data.size() - 1));
}

double calculateCV2(const std::vector<double>& t) {
    double sum_val = 0.0;
    for (int n = 0; n < t.size() - 1; n++) {
        double abs_err = abs(t[n + 1] - t[n]);
        double tot_interval = t[n + 1] + t[n];
        sum_val += 2.0 * (abs_err / tot_interval);
    }
    return sum_val / (t.size() - 1);
}

template<size_t N>
inline void single_neuron_rhs(double t, const State<N>& y, State<N>& dydt, const params& p, double V_ext_val) {
    const double x = y[0];
    const double T_prime = y[1];
    const double V = y[2];
    const double V_prime = y[3];

    const double dUdx = -((x - 0.1) * (x - 0.1) + 0.1) - 180 * exp(-((x - 0.8) * (x - 0.8)) / 0.01) + 0.2 * sqrt(10.0) * p.i0 - 100 * pow(x, 99); // Note: using x^17 for stability/speed

    const double inv_012 = 1.0 / 0.12;
    const double arg = std::clamp(x * inv_012, -100.0, 100.0);
    const double exp_pos = std::exp(arg);
    const double exp_neg = 1.0 / exp_pos;

    const double r = 0.5 * (exp_pos + exp_neg);
    const double r_inv = 1.0 / r;
    const double r_inv2 = r_inv * r_inv;
    const double r_prime = 0.5 * (exp_pos - exp_neg) * inv_012;

    dydt[0] = dUdx + 0.63 * V - p.qT * T_prime;
    dydt[1] = p.cT * ((2.0 * V * V_prime * r - V * V * r_prime) * r_inv2) - p.kappa * T_prime;
    dydt[2] = ((V_ext_val - (1.0 + (p.R_ext * r_inv)) * V)) / p.tc;
    dydt[3] = (((V * p.R_ext * r_prime) * r_inv2) - (1 + (p.R_ext * r_inv)) * V_prime) / p.tc;
}

template<size_t N, typename RHS>
inline void RK4_step(double& t, State<N>& y, const params& p, double dt, RHS&& rhs, double V_ext_val) {
    State<N> k1, k2, k3, k4, tmp;

    rhs(t, y, k1, p, V_ext_val);
    for (size_t i = 0; i < N; ++i) tmp[i] = y[i] + dt * 0.5 * k1[i];

    rhs(t + 0.5 * dt, tmp, k2, p, V_ext_val);
    for (size_t i = 0; i < N; ++i) tmp[i] = y[i] + dt * 0.5 * k2[i];

    rhs(t + 0.5 * dt, tmp, k3, p, V_ext_val);
    for (size_t i = 0; i < N; ++i) tmp[i] = y[i] + dt * k3[i];

    rhs(t + dt, tmp, k4, p, V_ext_val);
    for (size_t i = 0; i < N; ++i) y[i] += (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}

class SelectivitySimulator {
public:
    params p_global; // read only global parameters
    bool verbose;

    SelectivitySimulator(params initial_p, bool v = false)
        : p_global(initial_p), verbose(v) {}

    // parallel Batch Runner
    std::vector<SimulationOutput> run_batch(const std::vector<SimulationInput>& inputs) {
        std::vector<SimulationOutput> results(inputs.size());
        int completed = 0;
        int total = static_cast<int>(inputs.size());

        #pragma omp parallel for
        for (int i = 0; i < total; i++) {
            // Updated call: removed simulation_progress reference
            results[i] = core_logic(inputs[i], p_global, verbose);

            #pragma omp atomic
            completed++;

            // update progress bar
            #pragma omp critical
            {
                std::cout << "\rBatch Progress: " << (100 * completed / total) << "% [" 
                          << completed << "/" << total << " simulations]" << std::flush;
            }
        }
        std::cout << std::endl;
        return results;
    }

private:
    static SimulationOutput core_logic(const SimulationInput& in, const params& global_params, bool verbose) {
        // local copy of params 
        params p = global_params;
        p.V_DC = in.V_DC;
        p.R_ext = in.R_ext;

        // initial conditions
        double x0 = -1 + 0.00001;
        double T_prime0 = 0.0;
        double V0 = 0.0;
        double V_prime0 = 0.0;
        State<D> initial_state = { x0, T_prime0, V0, V_prime0 };

        double f_min = std::clamp(p.target_freq - p.frequency_range, 0.001, 1e5);
        double f_max = std::clamp(p.target_freq + p.frequency_range, 0.001, 1e5);
        int f_steps = static_cast<int>(std::round(std::abs(f_max - f_min) / std::abs(p.frequency_step)));

        std::vector<ParametricSweepResult> sweep_results(f_steps + 1);

        if (verbose) {
#pragma omp critical
            {
                std::cout << "[Thread " << omp_get_thread_num()
                    << "] Starting Sim: V_DC=" << p.V_DC
                    << ", R_ext=" << p.R_ext << std::endl;
            }
        }
        for (int i = 0; i <= f_steps; ++i) {
            State<D> current_state = initial_state;

            bool have_we_spiked = false;
            int spike_count = 0;
            double t = 0.0;
            int t_counter = 0;

            std::vector<double> spike_timings;
            spike_timings.reserve(static_cast<int>(p.max_spikes*2));

            double driving_frequency = f_min + (i * p.frequency_step);
            double angular_frequency = 2 * M_PI * driving_frequency;

            double V_ext_val = 0.0;

            // warmup
            while (spike_count < p.warmUp_spikes) {

                t = t_counter * p.dt;
                V_ext_val = p.V_DC + p.V_AC * cos(angular_frequency * t);
                RK4_step<D>(t, current_state, p, p.dt, single_neuron_rhs<D>, V_ext_val);

                if (current_state[0] > p.spike_threshold && !have_we_spiked) {
                    spike_count++;
                    have_we_spiked = true;
                }
                else if (current_state[0] < p.spike_threshold && have_we_spiked) {
                    have_we_spiked = false;
                }
                t_counter++;
            }
            spike_count = 0;
            // measurement
            while (spike_count < p.max_spikes) {
                t = t_counter * p.dt;
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
                t_counter++;
            }

            ParametricSweepResult sweep_result;
            double isi_mean = 0.0;
            double isi_standard_deviation = 0.0;
            double CV1 = 0.0;
            double CV2 = 0.0;

            if (spike_timings.size() >= 2) {
                std::vector<double> isi(spike_timings.size() - 1);
                for (size_t j = 0; j < isi.size(); ++j) {
                    isi[j] = spike_timings[j + 1] - spike_timings[j];
                }
                isi_standard_deviation = calculateSD(isi, isi_mean);
                if (isi_mean > 0.0) {
                    CV1 = isi_standard_deviation / isi_mean;
                    CV2 = calculateCV2(isi);
                }
            }

            sweep_result.tmax = spike_timings[spike_count-1]-spike_timings[0];
            sweep_result.spike_count = spike_count;
            sweep_result.winding_number = 1 / (isi_mean * driving_frequency);
            sweep_result.driving_frequency = driving_frequency;
            sweep_result.CV1 = CV1;
            sweep_result.CV2 = CV2;
            sweep_result.firing_rate = spike_count / sweep_result.tmax;
            sweep_result.isi_mean = isi_mean;

            sweep_results[i] = sweep_result;
        }
           
        if (verbose) {
#pragma omp critical
            {
                std::cout << "\n--- Thread " << omp_get_thread_num()
                    << " Results (V_DC: " << p.V_DC << ", R_ext: " << p.R_ext << ") ---" << std::endl;
                std::cout << "driving_freq,isi_mean,firing_rate,winding_number,CV1,CV2" << std::endl;

                for (const auto& res : sweep_results) {
                    std::cout << res.driving_frequency << ","
                        << res.isi_mean << ","
                        << res.firing_rate << ","   
                        << res.winding_number << ","
                        << res.CV1 << ","
                        << res.CV2 << std::endl;    
                }
                std::cout << "---------------------------------------------------------" << std::endl;
            }
        }

        // metrics Extraction
        double width = 0.0;
        std::vector<double> locking_frequencies;
        std::vector<int> locking_indices;

        for (int i = 0; i <= f_steps; i++) {
            if (abs(sweep_results[i].winding_number - 1.0) < 0.05) {
                width += p.frequency_step;
                locking_frequencies.push_back(sweep_results[i].driving_frequency);
                locking_indices.push_back(i);
            }
        }

        SimulationOutput out;
        if (!locking_frequencies.empty()) {
            out.is_phase_lock = 1;
            out.centre_of_phase_locking = (locking_frequencies.front() + locking_frequencies.back()) / 2.0;

            int max_index = static_cast<int>(sweep_results.size()) - 1;

            int left_index_offset = static_cast<int>(std::round(-p.depth_checking_range / p.frequency_step));
            int left_depth_index = std::clamp(locking_indices.front() + left_index_offset, 0, max_index);
            out.left_depth = sweep_results[left_depth_index].CV1;

            int right_index_offset = static_cast<int>(std::round(p.depth_checking_range / p.frequency_step));
            int right_depth_index = std::clamp(locking_indices.back() + right_index_offset, 0, max_index);
            out.right_depth = sweep_results[right_depth_index].CV1;

            out.depth = (out.left_depth + out.right_depth) / 2.0;
        }
        else {
            out.is_phase_lock = 0;
            out.width = 999.0;
            out.centre_of_phase_locking = 999.0;
            out.left_depth = 999.0;
            out.right_depth = 999.0;
            out.depth = 999.0;
        }

        out.frequency = 1/(sweep_results[0].isi_mean);
        out.width = width;
        out.CV1_at_target_freq = sweep_results[f_steps / 2].CV1;
        out.CV2_at_target_freq = sweep_results[f_steps / 2].CV2;

        return out;
    }
};

// debugging
int main() {
    std::cout << "Initializing Simulator..." << std::endl;

    params default_params;
    default_params.max_spikes = 10.0;
    default_params.warmUp_spikes = 10.0;
    bool verbose = true;
    SelectivitySimulator sim(default_params, verbose);

    // dummy batch 
    std::vector<SimulationInput> batch = {
        {6.0, 500.0},
        {6.1, 500.0},
        {6.2, 500.0},
        {6.3, 500.0},
        {6.4, 500.0},
        {6.5, 500.0},
        {6.5, 600.0},
        {7.0, 700.0},
        {7.5, 800.0}
    };

    std::cout << "Running batch of " << batch.size() << " parameter sets in parallel..." << std::endl;

    std::vector<SimulationOutput> results = sim.run_batch(batch);

    std::cout << "Simulation Finished! Results:" << std::endl;
    for (size_t i = 0; i < results.size(); ++i) {
        std::cout << "Run " << i << " | Locked: " << results[i].is_phase_lock
            << " | Width: " << results[i].width
            << " | Depth: " << results[i].depth << std::endl;
    }

    return 0;
}


// Instructiong for pybind11 linking to python
PYBIND11_MODULE(SelectivityLib, m) {
    m.doc() = "Parallel Neuron Selectivity Simulator";

    // params struct 
    py::class_<params>(m, "Params")
        .def(py::init<>())
        .def_readwrite("V_DC", &params::V_DC)
        .def_readwrite("R_ext", &params::R_ext)
        .def_readwrite("target_freq", &params::target_freq)
        .def_readwrite("V_AC", &params::V_AC)
        .def_readwrite("max_spikes", &params::max_spikes)
        .def_readwrite("warmUp_spikes", &params::warmUp_spikes)
        .def_readwrite("dt", &params::dt)
        .def_readwrite("frequency_range", &params::frequency_range)
        .def_readwrite("frequency_step", &params::frequency_step);

    // input struct
    py::class_<SimulationInput>(m, "SimulationInput")
        .def(py::init<double, double>())
        .def_readwrite("V_DC", &SimulationInput::V_DC)
        .def_readwrite("R_ext", &SimulationInput::R_ext);

    // output struct
    py::class_<SimulationOutput>(m, "SimulationOutput")
        .def_readonly("is_phase_lock", &SimulationOutput::is_phase_lock)
        .def_readonly("frequency", &SimulationOutput::frequency)
        .def_readonly("width", &SimulationOutput::width)
        .def_readonly("CV1_at_target_freq", &SimulationOutput::CV1_at_target_freq)
        .def_readonly("CV2_at_target_freq", &SimulationOutput::CV2_at_target_freq)
        .def_readonly("centre_of_phase_locking", &SimulationOutput::centre_of_phase_locking)
        .def_readonly("left_depth", &SimulationOutput::left_depth)
        .def_readonly("right_depth", &SimulationOutput::right_depth)
        .def_readonly("depth", &SimulationOutput::depth);

    py::class_<SelectivitySimulator>(m, "Simulator")
        .def(py::init<params, bool>(), py::arg("initial_p"), py::arg("verbose")=false)
        .def_readwrite("verbose", &SelectivitySimulator::verbose)
        .def("run_batch", &SelectivitySimulator::run_batch, "Runs a batch of inputs in parallel");
}
