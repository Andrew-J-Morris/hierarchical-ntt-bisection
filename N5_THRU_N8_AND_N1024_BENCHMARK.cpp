#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace std::chrono;

// ============================================================================
// PART 1: HARDWARE-NATIVE ENGINES FOR N = 5, 6, 7, 8 (TABLE II BENCHMARKS)
// ============================================================================
void evaluate_bisection_benchmarks(long long R) {
    long long R2 = R * R;
    auto t0 = high_resolution_clock::now();

    vector<long long> f1(R2 + 1, 0);
    f1[0] = 1;
    for (long long k = 1; k <= R; ++k) {
        if (k * k <= R2) f1[k * k] = 2;
    }

    vector<long long> f2(R2 + 1, 0);
    for (long long x = -R; x <= R; ++x) {
        long long x2 = x * x;
        for (long long y = -R; y <= R; ++y) {
            long long s = x2 + y * y;
            if (s <= R2) f2[s]++;
        }
    }

    vector<long long> f4(R2 + 1, 0);
    for (long long u = 0; u <= R2; ++u) {
        if (f2[u] == 0) continue;
        for (long long v = 0; v <= R2 - u; ++v) {
            if (f2[v] != 0) f4[u + v] += f2[u] * f2[v];
        }
    }

    vector<long long> F4(R2 + 1, 0);
    long long run_sum = 0;
    for (long long m = 0; m <= R2; ++m) {
        run_sum += f4[m];
        F4[m] = run_sum;
    }

    long long N5 = F4[R2];
    for (long long x5 = 1; x5 <= R; ++x5) {
        long long rem = R2 - x5 * x5;
        if (rem >= 0) N5 += 2 * F4[rem];
    }

    long long N6 = 0;
    for (long long v = 0; v <= R2; ++v) {
        if (f2[v] != 0) N6 += f2[v] * F4[R2 - v];
    }

    vector<long long> f6(R2 + 1, 0);
    for (long long u = 0; u <= R2; ++u) {
        if (f2[u] == 0) continue;
        for (long long v = 0; v <= R2 - u; ++v) {
            if (f4[v] != 0) f6[u + v] += f2[u] * f4[v];
        }
    }

    vector<long long> F6(R2 + 1, 0);
    run_sum = 0;
    for (long long m = 0; m <= R2; ++m) {
        run_sum += f6[m];
        F6[m] = run_sum;
    }

    long long N7 = F6[R2];
    for (long long x7 = 1; x7 <= R; ++x7) {
        long long rem = R2 - x7 * x7;
        if (rem >= 0) N7 += 2 * F6[rem];
    }

    unsigned long long N8 = 0;
    for (long long v = 0; v <= R2; ++v) {
        if (f4[v] != 0) N8 += f4[v] * F4[R2 - v];
    }

    auto t1 = high_resolution_clock::now();
    double elapsed_ms = duration_cast<microseconds>(t1 - t0).count() / 1000.0;

    cout << "Radius R = " << setw(3) << R << " | Time: " << fixed << setprecision(3) << setw(8) << elapsed_ms << " ms" << endl;
    cout << "  N_5(" << setw(3) << R << ") = " << N5 << endl;
    cout << "  N_6(" << setw(3) << R << ") = " << N6 << endl;
    cout << "  N_7(" << setw(3) << R << ") = " << N7 << endl;
    cout << "  N_8(" << setw(3) << R << ") = " << N8 << endl;
    cout << "----------------------------------------------------------------------" << endl;
}

// ============================================================================
// PART 2: 1024-DIMENSIONAL NTT BISECTION HARDWARE ENGINE (N = 1024, R = 100)
// ============================================================================
const long long MOD = 998244353;
const long long G = 3;

long long modPow(long long base, long long exp) {
    long long res = 1;
    base %= MOD;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % MOD;
        base = (base * base) % MOD;
        exp /= 2;
    }
    return res;
}

long long modInverse(long long n) {
    return modPow(n, MOD - 2);
}

void ntt_transform(vector<long long>& a, bool invert) {
    int n = static_cast<int>(a.size());
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }
    for (int len = 2; len <= n; len <<= 1) {
        long long wlen = modPow(G, (MOD - 1) / len);
        if (invert) wlen = modInverse(wlen);
        for (int i = 0; i < n; i += len) {
            long long w = 1;
            for (int j = 0; j < len / 2; j++) {
                long long u = a[i + j];
                long long v = (a[i + j + len / 2] * w) % MOD;
                a[i + j] = (u + v < MOD ? u + v : u + v - MOD);
                a[i + j + len / 2] = (u - v >= 0 ? u - v : u - v + MOD);
                w = (w * wlen) % MOD;
            }
        }
    }
    if (invert) {
        long long n_inv = modInverse(n);
        for (long long& x : a) x = (x * n_inv) % MOD;
    }
}

void run_1024d_benchmark() {
    long long R = 100;
    long long R2 = R * R;
    int M = 16384;

    cout << "\n======================================================================" << endl;
    cout << "  PART 2: 1024D HYPERBALL NTT BISECTION HARDWARE BENCHMARK           " << endl;
    cout << "======================================================================" << endl;

    auto t0 = high_resolution_clock::now();

    vector<long long> A(M, 0);
    A[0] = 1;
    for (long long k = 1; k <= R; ++k) {
        if (k * k <= R2) A[k * k] = 2;
    }

    ntt_transform(A, false);
    for (int i = 0; i < M; ++i) {
        A[i] = modPow(A[i], 1024);
    }
    ntt_transform(A, true);

    long long total_mod = 0;
    for (int m = 0; m <= R2; ++m) {
        total_mod = (total_mod + A[m]) % MOD;
    }

    auto t1 = high_resolution_clock::now();
    double elapsed_ms = duration_cast<microseconds>(t1 - t0).count() / 1000.0;

    cout << "1024D Lattice Count (mod 998244353) : " << total_mod << endl;
    cout << "Hardware Benchmark Execution Time   : " << fixed << setprecision(3) << elapsed_ms << " ms" << endl;
    cout << "----------------------------------------------------------------------" << endl;
}

int main() {
    cout << "================================================----------------------" << endl;
    cout << "  HIERARCHICAL BISECTION HARDWARE ENGINE: BENCHMARK AUDIT             " << endl;
    cout << "================================================----------------------" << endl;

    vector<long long> test_radii = { 10, 45, 100, 200 };
    for (long long R : test_radii) {
        evaluate_bisection_benchmarks(R);
    }

    run_1024d_benchmark();
    return 0;
}