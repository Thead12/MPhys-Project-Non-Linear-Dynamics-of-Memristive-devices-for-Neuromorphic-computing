// rossler_pseudo_arclength.cpp
#define AUTODIFF_EIGEN_FOUND

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <limits>
#include <array>
#include <vector>
#include <cmath>

#include <Eigen/Dense>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

using namespace std;
using namespace autodiff;

static struct params {
    double qT = 1.0;
    double kappa = 0.9;
    double R_ext = 500;
    double cT = 0.18;
    double i0 = 0.0;

    double lam = 0.12;
    double tc = 30.0;

    double V_ext = 0.0;
    double V_min = 0.0;
    double V_max = 0.0;
    
    double ds = 0.1;
    double newton_tol = 1e-8;
    int newton_maxiter = 20;
};

constexpr int D = 4;

using Scalar = dual;
template<int N>
using State = Eigen::Matrix<Scalar, N, 1>;

template<int N>
using StateDouble = Eigen::Matrix<double, N, 1>;

template<int N, typename T>
Eigen::Matrix<T, N, 1> single_neuron_rhs_templ(const Eigen::Matrix<T, N, 1>& y, const params& p, const T& V_ext);

template<int N>
void build_extended_z(const Eigen::Matrix<double, N, 1>& x_d, double Vext, Eigen::Matrix<Scalar, N + 1, 1>& z_ad);

template<int N>
bool newton_solve_equilibrium(StateDouble<N>& x0, params& p);

template<int N>
void compute_Jx_and_Jp(const StateDouble<N>& x_d, const params& p,
    Eigen::Matrix<double, N, N>& Jx_out, Eigen::Matrix<double, N, 1>& Jp_out,
    Eigen::Matrix<double, N, 1>& F_out);

template<int N>
void tangent_prediction_from_J(const Eigen::Matrix<double, N, N>& J_x,
    const Eigen::Matrix<double, N, 1>& J_p,
    Eigen::Matrix<double, N + 1, 1>& tangent_vec);

template<int N>
bool pseudo_arclength_step(Eigen::Matrix<double, N, 1>& x_curr,
    double& lambda_curr,
    const Eigen::Matrix<double, N + 1, 1>& tangent_curr,
    params& p,
    Eigen::Matrix<double, N + 1, 1>& tangent_next,
    StateDouble<N>& x_next_out,
    double& lambda_next_out);

template<int N, typename T>
Eigen::Matrix<T, N, 1> single_neuron_rhs_templ(const Eigen::Matrix<T, N, 1>& y, const params& p, const T& V_ext)
{
    const T x = y(0);
    const T T_prime = y(1);
    const T V = y(2);
    const T V_prime = y(3);

    const T dUdx = -((x - T(0.1)) * (x - T(0.1)) + T(0.1))
        - T(180) * exp(-((x - T(0.8)) * (x - T(0.8))) / T(0.01))
        + T(0.2) * sqrt(T(10.0)) * T(p.i0)
        - T(100.0) * pow(x, 99);

    const T r = cosh(x / T(0.12));
    const T r2 = r * r;
    const T r_prime = sinh(x / T(0.12)) / T(0.12);

    Eigen::Matrix<T, N, 1> dydt;
    dydt(0) = dUdx + T(0.63) * V - T(p.qT) * T_prime;
    dydt(1) = T(p.cT) * ((T(2.0) * V * V_prime * r - V * V * r_prime) / r2) - T(p.kappa) * T_prime;
    dydt(2) = ((V_ext - (T(1.0) + T(p.R_ext) / r) * V)) / T(p.tc);
    dydt(3) = ((V * T(p.R_ext) * r_prime) / r2 - (T(1.0) + T(p.R_ext) / r) * V_prime) / T(p.tc);

    return dydt;
}

// ------------------------------- helpers -------------------------------

template<int N>
void build_extended_z(const Eigen::Matrix<double, N, 1>& x_d, double Vext, Eigen::Matrix<Scalar, N + 1, 1>& z_ad)
{
    for (int i = 0; i < N; i++) z_ad(i) = Scalar(x_d(i));
    z_ad(N) = Scalar(Vext); // last element is lambda
}

template<int N>
bool newton_solve_equilibrium(StateDouble<N>& x0, params& p)
{
    for (int iter = 0; iter < p.newton_maxiter; ++iter)
    {
        Eigen::Matrix<double, N, N> Jx;
        Eigen::Matrix<double, N, 1> Jp_dummy;
        Eigen::Matrix<double, N, 1> Fd;
        compute_Jx_and_Jp<N>(x0, p, Jx, Jp_dummy, Fd); // Fd = F(x0,p.V_ext)

        double resnorm = Fd.template lpNorm<2>(); // norm of Fd
        if (resnorm < p.newton_tol) return true;

        // solve Jx * dx = -F
        Eigen::PartialPivLU<Eigen::Matrix<double, N, N>> solver(Jx);
        Eigen::Matrix<double, N, 1> dx = solver.solve(-Fd);

        x0 += dx;

        if (dx.template lpNorm<2>() < p.newton_tol) return true;
    }
    return false;
}

template<int N>
void compute_Jx_and_Jp(const StateDouble<N>& x_d, const params& p,
    Eigen::Matrix<double, N, N>& Jx_out, Eigen::Matrix<double, N, 1>& Jp_out,
    Eigen::Matrix<double, N, 1>& F_out)
{
    // compute J_x and F
    Eigen::Matrix<Scalar, N, 1> x_ad;
    for (int i = 0; i < N; ++i) x_ad(i) = Scalar(x_d(i));

    auto f_wrt_x = [&](const Eigen::Matrix<Scalar, N, 1>& xa)->Eigen::Matrix<Scalar, N, 1>
    {
        return single_neuron_rhs_templ<N, Scalar>(xa, p, Scalar(p.V_ext));
    };

    Eigen::Matrix<Scalar, N, 1> F_ad;
    Eigen::Matrix<Scalar, N, N> Jx_ad;
    Jx_ad = jacobian(f_wrt_x, wrt(x_ad), at(x_ad), F_ad);

    // convert to doubles
    for (int i = 0; i < N; ++i) {
        F_out(i) = val(F_ad(i));
        for (int j = 0; j < N; ++j) {
            Jx_out(i, j) = val(Jx_ad(i, j));
        }
    }

    // compute J_p
    Scalar Vext_ad = Scalar(p.V_ext);
    auto f_wrt_V = [&](const Scalar& Vext)->Eigen::Matrix<Scalar, N, 1>
    {
        return single_neuron_rhs_templ<N, Scalar>(x_ad, p, Vext);
    };

    Eigen::Matrix<Scalar, N, 1> Jp_ad;
    Eigen::Matrix<Scalar, N, 1> Fp_ad; // unused, just there for jacobian function paramters
    Jp_ad = jacobian(f_wrt_V, wrt(Vext_ad), at(Vext_ad), Fp_ad);

    for (int i = 0; i < N; i++) Jp_out(i) = val(Jp_ad(i)); // converting to double
}

// scd decomposition method
template<int N>
void tangent_prediction_from_J(const Eigen::Matrix<double, N, N>& J_x,
    const Eigen::Matrix<double, N, 1>& J_p,
    Eigen::Matrix<double, N + 1, 1>& tangent_vec)
{
    Eigen::Matrix<double, N, N + 1> A;
    for (int i = 0; i < N; i++) {
        for (int j=0;j<N;j++) {
            A(i, j) = J_x(i, j);
        }
    }
    for (int k = 0; k < N; k++) {
        A(k, N) = J_p(k);
    }

    Eigen::JacobiSVD < Eigen::Matrix<double, N, N + 1>> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Matrix<double, N+1,N+1> V = svd.matrixV();
    Eigen::Matrix<double, N + 1, 1> v = V.col(V.cols() - 1);

    cerr << "\nt_x(0) normalised = " << tangent_vec(0) << endl;
    cerr << "t_lam normalised = " << tangent_vec(N) << endl;

    if (v.dot(tangent_vec) < 0.0) {
        tangent_vec = -v.normalized();
    }
    tangent_vec = v.normalized();
}

// old method
//template<int N>
//void tangent_prediction_from_J(const Eigen::Matrix<double, N, N>& J_x,
//    const Eigen::Matrix<double, N, 1>& J_p,
//    Eigen::Matrix<double, N + 1, 1>& tangent_vec)
//{
//    double t_lam = 1.0;
//    Eigen::Matrix<double, N, 1> t_x = -J_x.colPivHouseholderQr().solve(J_p * t_lam);
//
//    for (int i = 0; i < N; ++i) tangent_vec(i) = t_x(i);
//    tangent_vec(N) = t_lam;
//    tangent_vec.normalize();
//
//    cerr << "\nt_x(0) normalised = " << tangent_vec(0) << endl;
//    cerr << "t_lam normalised = " << tangent_vec(N) << endl;
//}

template<int N>
bool pseudo_arclength_step(
    Eigen::Matrix<double, N, 1>& x_curr,
    double& lambda_curr,
    const Eigen::Matrix<double, N + 1, 1>& tangent_curr,
    params& p,
    Eigen::Matrix<double, N + 1, 1>& tangent_next,
    StateDouble<N>& x_next_out,
    double& lambda_next_out)
{
    // Predictor: X_pred = X_curr + ds * tangent_curr
    Eigen::Matrix<double, N + 1, 1> Xcurr_ext;
    for (int i = 0; i < N; ++i) Xcurr_ext(i) = x_curr(i);
    Xcurr_ext(N) = lambda_curr;

    Eigen::Matrix<double, N + 1, 1> Xpred_ext = Xcurr_ext + p.ds * tangent_curr;

    Eigen::Matrix<double, N, 1> x_pred;
    double lambda_pred;
    for (int i = 0; i < N; i++) x_pred(i) = Xpred_ext(i);
    lambda_pred = Xpred_ext(N);

    // Corrector: H(X) = [ F(x,λ); t^T (X - X_curr) - ds ] = 0
    Eigen::Matrix<double, N + 1, 1> X_ext;
    for (int i = 0; i < N; ++i) X_ext(i) = x_pred(i);
    X_ext(N) = lambda_pred;

    // Newton method
    for (int iter = 0; iter < p.newton_maxiter; ++iter)
    {
        // acobians
        Eigen::Matrix<double, N, N> Jx;
        Eigen::Matrix<double, N, 1> Jp;
        Eigen::Matrix<double, N, 1> Fd;

        params ptemp = p;
        ptemp.V_ext = X_ext(N);

        // compute Jx, Jp, F at (X_ext(0:N-1), ptemp.V_ext)
        StateDouble<N> x_temp;
        for (int i = 0; i < N; ++i) x_temp(i) = X_ext(i);

        compute_Jx_and_Jp<N>(x_temp, ptemp, Jx, Jp, Fd);

        // H = [F; constraint]
        Eigen::Matrix<double, N + 1, 1> H;
        for (int i = 0; i < N; ++i) H(i) = Fd(i);

        // constraint c = tangent_curr^T (X_ext - Xcurr_ext) - ds
        double c = 0.0;
        for (int i = 0; i < N + 1; i++) c += tangent_curr(i) * (X_ext(i) - Xcurr_ext(i));
        c -= p.ds;
        H(N) = c;

        double resnorm = H.template lpNorm<2>();
        if (resnorm < p.newton_tol) // if converged
        {
            for (int i = 0; i < N; i++) x_next_out(i) = X_ext(i);
            lambda_next_out = X_ext(N);

            // compute Jacobians at new point to get tangent_next
            Eigen::Matrix<double, N, N> Jx_new;
            Eigen::Matrix<double, N, 1> Jp_new;
            Eigen::Matrix<double, N, 1> Fnew;
            params pnew = p;
            pnew.V_ext = lambda_next_out;
            compute_Jx_and_Jp<N>(x_next_out, pnew, Jx_new, Jp_new, Fnew);

            tangent_prediction_from_J<N>(Jx_new, Jp_new, tangent_next);
            return true;
        }

        // augmented jacobian
        Eigen::Matrix<double, N + 1, N + 1> Jaug;
        Jaug.setZero();
        Jaug.template block<N, N>(0, 0) = Jx;
        Jaug.template block<N, 1>(0, N) = Jp;
        for (int j = 0; j < N; ++j) Jaug(N, j) = tangent_curr(j);
        Jaug(N, N) = tangent_curr(N);

        // solve Jaug * dX = -H
        Eigen::PartialPivLU<Eigen::Matrix<double, N + 1, N + 1>> solver(Jaug);
        Eigen::Matrix<double, N + 1, 1> dX = solver.solve(-H);

        for (int i = 0; i < N + 1; ++i) X_ext(i) += dX(i);
    }

    // if reached here, Newton did not converge
    return false;
}

// ============================== main ==============================
int main(int argc, char* argv[])
{
    if (argc < 3) {
        cerr << "Usage: single_neuron.exe <params_file> <output_csv>\n";
        return 1;
    }
    string param_file = argv[1];
    string output_file = argv[2];

    ifstream infile(param_file);
    if (!infile.is_open()) {
        cerr << "Failed to open " << param_file << "\n";
        return 1;
    }

    ofstream outfile(output_file);
    if (!outfile.is_open()) {
        cerr << "Failed to open " << output_file << " for writing\n";
        return 1;
    }

    // initial conditions
    double x0 = -1 + 0.00001;
    double T_prime0 = 0.0;
    double V0 = 0.0;
    double V_prime0 = 0.0;
    StateDouble<D> state_d;
    state_d << x0, T_prime0, V0, V_prime0;

    params p;
    double tmax = 0.0, dt = 0.0;

    // read params file
    string line;
    while (getline(infile, line)) {
        auto comment_pos = line.find('#');
        if (comment_pos != string::npos) line = line.substr(0, comment_pos);
        if (line.empty()) continue;
        istringstream iss(line);
        string name; double value;
        if (!(iss >> name >> value)) continue;
        if (name == "x0") state_d(0) = value;
        else if (name == "T_prime0") state_d(1) = value;
        else if (name == "V0") state_d(2) = value;
        else if (name == "V_prime0") state_d(3) = value;
        else if (name == "qT") p.qT = value;
        else if (name == "kappa") p.kappa = value;
        else if (name == "R_ext") p.R_ext = value;
        else if (name == "cT") p.cT = value;
        else if (name == "i0") p.i0 = value;
        else if (name == "lam") p.lam = value;
        else if (name == "tc") p.tc = value;
        else if (name == "V_min") p.V_min = value;
        else if (name == "V_max") p.V_max = value;
        else if (name == "ds") p.ds = value;
        else if (name == "newton_tol") p.newton_tol = value;
        else if (name == "newton_maxiter") p.newton_maxiter = value;
        else if (name == "tmax") tmax = value;
        else if (name == "dt") dt = value;

        else cerr << "Unknown parameter: " << name << "\n";
    }

    p.V_ext = p.V_min;

    // write header
    outfile << "Parameters: V_min=" << p.V_min << ", V_max=" << p.V_max << ", qT=" << p.qT << ", kappa=" << p.kappa
        << ", R_ext=" << p.R_ext << ", cT=" << p.cT << ", i0=" << p.i0
        << ", tc=" << p.tc << ", lam=" << p.lam
        << ", tmax=" << tmax << ", dt=" << dt << endl;
    outfile << "V_ext";
    for (int i = 0; i < D; ++i) outfile << ", x" << i;
    outfile << endl;

    // initial equilibrium at V_ext = V_min using Newton
    StateDouble<D> x_eq = state_d;
    bool ok = newton_solve_equilibrium<D>(x_eq, p);
    if (!ok) {
        cerr << "Initial Newton to find equilibrium failed at V_ext = " << p.V_ext << endl;
        return 1;
    }

    // initial Jacobians and tangent
    Eigen::Matrix<double, D, D> Jx0;
    Eigen::Matrix<double, D, 1> Jp0;
    Eigen::Matrix<double, D, 1> F0;
    compute_Jx_and_Jp<D>(x_eq, p, Jx0, Jp0, F0);

    Eigen::Matrix<double, D + 1, 1> tangent0;
    tangent_prediction_from_J<D>(Jx0, Jp0, tangent0);

    outfile << p.V_ext;
    for (int i = 0; i < D; ++i) outfile << ", " << setprecision(15) << x_eq(i);
    outfile << endl;

    // continuation loop
    Eigen::Matrix<double, D + 1, 1> tangent_curr = tangent0;
    Eigen::Matrix<double, D, 1> x_curr = x_eq;
    double lambda_curr = p.V_ext;

    int max_steps = int(ceil((p.V_max - p.V_min) / fabs(p.ds))) + 1;
    for (int step = 0; step < max_steps; ++step)
    {
        // stopping condition
        if ((p.ds > 0 && lambda_curr >= p.V_max) || (p.ds < 0 && lambda_curr <= p.V_max)) break;

        // continuation step (predictor + corrector)
        StateDouble<D> x_next;
        double lambda_next;
        Eigen::Matrix<double, D + 1, 1> tangent_next;

        bool succ = false;
        double tmp_ds = p.ds;
        while (true)
        {
            cerr << "We in here" << endl;
            succ = pseudo_arclength_step<D>(x_curr, lambda_curr, tangent_curr, p, tangent_next, x_next, lambda_next);
            if (succ) {
                break;
            }
            cerr << "Continuation step failed at lambda = " << lambda_curr << ". " << endl;
            p.ds = p.ds / 10;
            cerr << "Reducing arclength ds = " << p.ds << endl;
        }
        p.ds = tmp_ds;
        

        x_curr = x_next;
        lambda_curr = lambda_next;
        tangent_curr = tangent_next;

        // print
        outfile << lambda_curr;
        for (int i = 0; i < D; ++i) outfile << ", " << setprecision(15) << x_curr(i);
        outfile << endl;

        cerr << "| step " << step << " lambda=" << lambda_curr << "\n";
    }

    cerr << "|=============== Continuation Finished ===============|\n";
    return 0;
}
