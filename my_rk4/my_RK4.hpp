double RK4(double t, double y0, double dt, double rhs(double t, double y));
vector<double> RK4(double t, vector<double> y0, double dt, vector<double> rhs(double t, vector<double> y, int n), int n);
vector<vector<double>> numerical_jacobian(double t, vector<double> x, double h, vector<double> rhs(double t, vector<double> x, int n), int n, int m);
double mse(vector<double> x, vector<double> y);
double system_mse(vector<vector<double>> x, vector<vector<double>> y, int n);