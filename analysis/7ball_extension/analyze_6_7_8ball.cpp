/**
 * Unified Hyperball Lattice Residual & Spectral Suite (Dimensions 6, 7, 8)
 * Targets     : b055412.txt (6D), b055413_NTT.txt (7D), b055414.txt (8D)
 * Environment : Microsoft Visual Studio (C++17 or later)
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <algorithm>
#include <numeric>

using namespace std;

constexpr double PI_VAL = 3.14159265358979323846;

// ============================================================================
// 1. EXACT CONTINUOUS VOLUME COEFFICIENTS
// ============================================================================
double get_volume_coefficient(int dimension) {
    switch (dimension) {
    case 6: return std::pow(PI_VAL, 3) / 6.0;             // pi^3 / 6
    case 7: return (16.0 * std::pow(PI_VAL, 3)) / 105.0;  // 16*pi^3 / 105
    case 8: return std::pow(PI_VAL, 4) / 24.0;            // pi^4 / 24
    default:
        cerr << "[-] Unsupported dimension: " << dimension << "\n";
        return 0.0;
    }
}

// ============================================================================
// 1. DATA INGESTION STRUCT (C26495 COMPLIANT)
// ============================================================================
struct LatticeData {
    int dimension = 0;            // In-class default initialization silences C26495
    vector<double> R;
    vector<double> N_R;
    vector<double> E_R;
};

LatticeData load_and_compute_residuals(const string& filepath, int dim) {
    LatticeData data;
    data.dimension = dim;
    double coeff = get_volume_coefficient(dim);

    ifstream file(filepath);
    if (!file.is_open()) {
        cerr << "[-] Error: Could not open file: " << filepath << "\n";
        return data;
    }

    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        double r_val;
        string count_str;
        if (ss >> r_val >> count_str) {
            if (r_val <= 0.0) continue;

            double n_val = stod(count_str);
            double v_val = coeff * std::pow(r_val, dim);
            double e_val = n_val - v_val;

            data.R.push_back(r_val);
            data.N_R.push_back(n_val);
            data.E_R.push_back(e_val);
        }
    }
    return data;
}

// ============================================================================
// 2. LOG-LOG LINEAR REGRESSION
// ============================================================================
struct RegressionResult {
    double alpha;
    double intercept;
    double r_squared;
};

RegressionResult log_log_regression(const vector<double>& R, const vector<double>& E, double r_min, double r_max) {
    vector<double> x, y;
    for (size_t i = 0; i < R.size(); ++i) {
        if (R[i] >= r_min && R[i] <= r_max && std::abs(E[i]) > 0.0) {
            x.push_back(std::log(R[i]));
            y.push_back(std::log(std::abs(E[i])));
        }
    }

    size_t n = x.size();
    if (n < 2) return { 0.0, 0.0, 0.0 };

    double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
    double sum_y = std::accumulate(y.begin(), y.end(), 0.0);
    double mean_x = sum_x / n;
    double mean_y = sum_y / n;

    double num = 0.0, den = 0.0;
    for (size_t i = 0; i < n; ++i) {
        num += (x[i] - mean_x) * (y[i] - mean_y);
        den += (x[i] - mean_x) * (x[i] - mean_x);
    }

    double slope = num / den;
    double intercept = mean_y - slope * mean_x;

    double ss_tot = 0.0, ss_res = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double y_pred = slope * x[i] + intercept;
        ss_res += (y[i] - y_pred) * (y[i] - y_pred);
        ss_tot += (y[i] - mean_y) * (y[i] - mean_y);
    }
    double r2 = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 1.0;

    return { slope, intercept, r2 };
}

// ============================================================================
// 3. QUADRATIC LEAST-SQUARES DETRENDING
// ============================================================================
vector<double> quadratic_detrend(const vector<double>& x, const vector<double>& y) {
    size_t n = x.size();
    double s0 = n, s1 = 0, s2 = 0, s3 = 0, s4 = 0;
    double t0 = 0, t1 = 0, t2 = 0;

    for (size_t i = 0; i < n; ++i) {
        double xi = x[i];
        double yi = y[i];
        double x2 = xi * xi;
        s1 += xi; s2 += x2; s3 += x2 * xi; s4 += x2 * x2;
        t0 += yi; t1 += yi * xi; t2 += yi * x2;
    }

    double A[3][4] = {
        {s4, s3, s2, t2},
        {s3, s2, s1, t1},
        {s2, s1, s0, t0}
    };

    for (int i = 0; i < 3; ++i) {
        int pivot = i;
        for (int j = i + 1; j < 3; ++j) {
            if (std::abs(A[j][i]) > std::abs(A[pivot][i])) pivot = j;
        }
        for (int k = 0; k < 4; ++k) std::swap(A[i][k], A[pivot][k]);

        double diag = A[i][i];
        for (int k = i; k < 4; ++k) A[i][k] /= diag;
        for (int j = 0; j < 3; ++j) {
            if (j != i) {
                double factor = A[j][i];
                for (int k = i; k < 4; ++k) A[j][k] -= factor * A[i][k];
            }
        }
    }

    double a = A[0][3], b = A[1][3], c = A[2][3];

    vector<double> detrended(n);
    for (size_t i = 0; i < n; ++i) {
        detrended[i] = y[i] - (a * x[i] * x[i] + b * x[i] + c);
    }
    return detrended;
}

// ============================================================================
// 4. COOLEY-TUKEY FFT
// ============================================================================
void fft(vector<complex<double>>& a) {
    size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    for (size_t len = 2; len <= n; len <<= 1) {
        double ang = -2.0 * PI_VAL / len;
        complex<double> wlen(std::cos(ang), std::sin(ang));
        for (size_t i = 0; i < n; i += len) {
            complex<double> w(1.0);
            for (size_t j = 0; j < len / 2; j++) {
                complex<double> u = a[i + j];
                complex<double> v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

// ============================================================================
// SUITE EXECUTION FUNCTION
// ============================================================================
void run_analysis(const string& filename, int dim) {
    cout << "\n======================================================================\n";
    cout << "    " << dim << "-BALL LATTICE RESIDUAL & SPECTRAL SUITE (R = 1..10000)\n";
    cout << "    Source File: " << filename << "\n";
    cout << "======================================================================\n";

    LatticeData data = load_and_compute_residuals(filename, dim);
    if (data.R.empty()) return;

    cout << "[+] Dataset Successfully Loaded: " << data.R.size() << " radial points.\n";
    cout << "[+] Maximum Radius: " << fixed << setprecision(0) << data.R.back() << "\n\n";

    // 1. Bracketed Log-Log Regression
    cout << "--- 1. Error Growth Exponent alpha (|E(R)| ~ C * R^alpha) ---\n";
    vector<pair<double, double>> brackets = {
        {100, 1000}, {1000, 2000}, {2000, 3000}, {3000, 4000}, {4000, 5000},
        {5000, 6000}, {6000, 7000}, {7000, 8000}, {8000, 9000}, {9000, 10000},
        {1000, 5000}, {5000, 10000}, {1000, 10000}
    };

    cout << left << setw(20) << "Radial Bracket"
        << setw(16) << "Exponent (alpha)"
        << setw(12) << "R^2 Fit" << "\n";
    cout << string(48, '-') << "\n";

    for (auto& b : brackets) {
        RegressionResult res = log_log_regression(data.R, data.E_R, b.first, b.second);
        stringstream ss_bracket;
        ss_bracket << "[" << (int)b.first << " - " << (int)b.second << "]";
        cout << left << setw(20) << ss_bracket.str()
            << fixed << setprecision(4) << setw(16) << res.alpha
            << setprecision(5) << setw(12) << res.r_squared << "\n";
    }

    // 2. High-End Spectral FFT Analysis (R in [1000, 10000])
    cout << "\n--- 2. Spectral FFT Modes on Detrended Residuals (R in [1000, 10000]) ---\n";
    vector<double> r_sub, e_sub;
    for (size_t i = 0; i < data.R.size(); ++i) {
        if (data.R[i] >= 1000.0) {
            r_sub.push_back(data.R[i]);
            e_sub.push_back(data.E_R[i]);
        }
    }

    vector<double> detrended = quadratic_detrend(r_sub, e_sub);
    size_t orig_len = detrended.size();
    size_t fft_len = 1;
    while (fft_len < orig_len) fft_len <<= 1;

    vector<complex<double>> signal(fft_len, 0.0);
    for (size_t i = 0; i < orig_len; ++i) signal[i] = detrended[i];
    fft(signal);

    size_t half_len = fft_len / 2;
    vector<pair<double, size_t>> spectrum;
    for (size_t i = 1; i < half_len; ++i) {
        spectrum.push_back({ std::norm(signal[i]), i });
    }
    std::sort(spectrum.rbegin(), spectrum.rend());

    cout << left << setw(8) << "Rank"
        << setw(18) << "Frequency (cyc/u)"
        << setw(20) << "Period (rad units)"
        << setw(16) << "Spectral Power" << "\n";
    cout << string(62, '-') << "\n";

    for (int k = 0; k < 4 && k < (int)spectrum.size(); ++k) {
        size_t idx = spectrum[k].second;
        double freq = (double)idx / (double)fft_len;
        double period = 1.0 / freq;
        cout << left << setw(8) << (k + 1)
            << fixed << setprecision(5) << setw(18) << freq
            << setprecision(2) << setw(20) << period
            << scientific << setprecision(2) << setw(16) << spectrum[k].first << "\n";
    }
    cout << "\n======================================================================\n";
}

int main() {
    // Sequentially processes 6D and 8D datasets
    run_analysis("b055412.txt", 6);
    run_analysis("b055414.txt", 8);
    return 0;
}