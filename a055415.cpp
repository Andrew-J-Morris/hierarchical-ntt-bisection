/**
 * 9-Ball NTT + CRT Exact Sequence Generator (OEIS A055415)
 * Framework   : 256-bit Mixed-Radix Accumulator (Plug-and-Play Integration)
 * Target      : R = 0 to 10000
 * Optimization: Full OpenMP Threading inside the NTT Butterfly Stages.
 * Memory Guard: Locked at ~8.58 GiB max footprint per active loop.
 */

#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <string>
#include <cstdint>
#include <algorithm>
#include <fstream>
#include <omp.h>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

using namespace std;
using namespace std::chrono;

// ============================================================================
// 1. HARDWARE-NATIVE ARITHMETIC INTRINSICS (MSVC x64 SAFE)
// ============================================================================
static inline uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t mod) {
#if defined(_MSC_VER) && defined(_M_X64)
    uint64_t hi, lo;
    lo = _umul128(a, b, &hi);
    uint64_t rem;
    _udiv128(hi, lo, mod, &rem);
    return rem;
#else
    return static_cast<uint64_t>((static_cast<unsigned __int128>(a) * b) % mod);
#endif
}

uint64_t modPow(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t res = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) res = mulmod64(res, base, mod);
        base = mulmod64(base, base, mod);
        exp >>= 1;
    }
    return res;
}

uint64_t modInverse(uint64_t n, uint64_t mod) {
    return modPow(n, mod - 2ULL, mod);
}

// ============================================================================
// 2. 256-BIT "GOLDILOCKS" ACCUMULATOR & STRING CONVERTER
// ============================================================================
struct uint256_t {
    uint64_t w[4] = { 0, 0, 0, 0 }; // lo, gl_lo, gl_hi, hi

    void mul_add_64(uint64_t multiplier, uint64_t addend) {
        uint64_t carry = addend;
        for (int i = 0; i < 4; ++i) {
            uint64_t high;
            uint64_t low = _umul128(w[i], multiplier, &high);
            unsigned char c1 = _addcarry_u64(0, low, carry, &w[i]);
            _addcarry_u64(c1, high, 0, &carry);
        }
    }

    uint32_t divmod32(uint32_t divisor) {
        uint64_t rem = 0;
        for (int i = 3; i >= 0; --i) {
            uint64_t cur = (rem << 32) | (w[i] >> 32);
            uint64_t q_high = cur / divisor;
            rem = cur % divisor;

            cur = (rem << 32) | (w[i] & 0xFFFFFFFFULL);
            uint64_t q_low = cur / divisor;
            rem = cur % divisor;

            w[i] = (q_high << 32) | q_low;
        }
        return static_cast<uint32_t>(rem);
    }

    bool is_zero() const {
        return (w[0] | w[1] | w[2] | w[3]) == 0;
    }

    string to_string() const {
        if (is_zero()) return "0";
        uint256_t temp = *this;
        string s;
        while (!temp.is_zero()) {
            s += std::to_string(temp.divmod32(10));
        }
        reverse(s.begin(), s.end());
        return s;
    }
};

// ============================================================================
// 3. DYNAMIC NTT PRIME GENERATOR 
// ============================================================================
void generate_ntt_primes(size_t M, int count, vector<uint64_t>& primes, vector<uint64_t>& roots) {
    uint64_t k_start = (1ULL << 45) / static_cast<uint64_t>(M);
    if (k_start % 2 == 0) k_start++;
    uint64_t k = k_start;

    while (static_cast<int>(primes.size()) < count) {
        uint64_t p = k * static_cast<uint64_t>(M) + 1ULL;
        bool is_prime = true;

        if (p % 3 == 0 || p % 5 == 0 || p % 7 == 0 || p % 11 == 0 || p % 13 == 0) {
            is_prime = false;
        }
        else {
            for (uint64_t i = 17; i * i <= p; i += 2) {
                if (p % i == 0) {
                    is_prime = false;
                    break;
                }
            }
        }

        if (is_prime) {
            for (uint64_t g = 2; g < p; ++g) {
                if (modPow(g, (p - 1ULL) / 2ULL, p) == p - 1ULL) {
                    primes.push_back(p);
                    roots.push_back(g);
                    break;
                }
            }
        }
        k += 2;
    }
}

// ============================================================================
// 4. FULLY THREADED IN-PLACE NTT & 9-BALL EXTRACTION 
// ============================================================================
void ntt_transform(vector<uint64_t>& a, bool invert, uint64_t mod, uint64_t g) {
    int64_t n = static_cast<int64_t>(a.size());

    for (int64_t i = 1, j = 0; i < n; i++) {
        int64_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }

    for (int64_t len = 2; len <= n; len <<= 1) {
        uint64_t wlen = modPow(g, (mod - 1ULL) / static_cast<uint64_t>(len), mod);
        if (invert) wlen = modInverse(wlen, mod);

#pragma omp parallel for schedule(dynamic, 64)
        for (int64_t i = 0; i < n; i += len) {
            uint64_t w = 1;
            int64_t half_len = len / 2;
            for (int64_t j = 0; j < half_len; j++) {
                uint64_t u = a[i + j];
                uint64_t v = mulmod64(a[i + j + half_len], w, mod);
                a[i + j] = (u + v < mod) ? (u + v) : (u + v - mod);
                a[i + j + half_len] = (u >= v) ? (u - v) : (u + mod - v);
                w = mulmod64(w, wlen, mod);
            }
        }
    }

    if (invert) {
        uint64_t n_inv = modInverse(static_cast<uint64_t>(n), mod);
#pragma omp parallel for schedule(static)
        for (int64_t i = 0; i < n; i++) {
            a[i] = mulmod64(a[i], n_inv, mod);
        }
    }
}

void compute_9ball_sequence(uint64_t R, size_t M, uint64_t mod, uint64_t g, vector<uint64_t>& seq_out) {
    uint64_t R2 = R * R;

    // Scoped Buffer allocation - Keeps RAM usage strictly bounded.
    vector<uint64_t> A(M, 0);

    A[0] = 1;
    for (uint64_t k = 1; k <= R; ++k) {
        A[k * k] = 2;
    }

    cout << "    -> Forward NTT... " << flush;
    ntt_transform(A, false, mod, g);

    cout << "Done.\n    -> Parallel Exponentiation to Power 9... " << flush;
#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < static_cast<int64_t>(M); ++i) {
        A[i] = modPow(A[i], 9ULL, mod);
    }

    cout << "Done.\n    -> Inverse NTT... " << flush;
    ntt_transform(A, true, mod, g);

    cout << "Done.\n    -> Extracting 9D Shell Totals... " << flush;
    uint64_t total_mod = 0;
    uint64_t target_r = 0;
    uint64_t next_stop = 0;

    for (size_t m = 0; m <= R2; ++m) {
        total_mod = (total_mod + A[m]);
        if (total_mod >= mod) total_mod -= mod;

        if (m == next_stop) {
            seq_out[target_r] = total_mod;
            target_r++;
            if (target_r <= R) {
                next_stop = target_r * target_r;
            }
        }
    }
    cout << "Done.\n" << flush;
}

// ============================================================================
// 5. 256-BIT GARNER'S ALGORITHM (MIXED-RADIX CRT)
// ============================================================================
string solve_crt_garner_256(const vector<uint64_t>& p, const vector<uint64_t>& r) {
    size_t k = p.size();
    vector<uint64_t> v(k);

    for (size_t i = 0; i < k; ++i) {
        uint64_t cur = r[i];
        for (size_t j = 0; j < i; ++j) {
            uint64_t p_mod = p[i];
            uint64_t vj_mod = v[j] % p_mod;
            uint64_t diff = (cur >= vj_mod) ? (cur - vj_mod) : (cur + p_mod - vj_mod);
            uint64_t inv = modInverse(p[j] % p_mod, p_mod);
            cur = mulmod64(diff, inv, p_mod);
        }
        v[i] = cur;
    }

    uint256_t result;
    result.w[0] = v[k - 1];

    for (int i = static_cast<int>(k) - 2; i >= 0; --i) {
        result.mul_add_64(p[i], v[i]);
    }

    return result.to_string();
}

// ============================================================================
// MAIN DRIVER & FILE OUTPUT
// ============================================================================
int main() {
    uint64_t R = 10000;
    uint64_t R2 = R * R;

    // Acyclic constraint for 9-Ball: M >= 9 * R^2 + 1
    size_t M = 1;
    while (M <= 9ULL * R2) {
        M <<= 1;
    }

    vector<uint64_t> primes;
    vector<uint64_t> roots;

    // Generate FIVE 45-bit primes (225 bits composite capacity)
    generate_ntt_primes(M, 5, primes, roots);

    cout << "======================================================================\n";
    cout << " 9-BALL NTT + CRT EXACT SEQUENCE GENERATOR (OEIS A055415)\n";
    cout << " Target Radius R = " << R << " | Dynamic Transform Size M = " << M << "\n";
    cout << " Max Active RAM  = " << (M * sizeof(uint64_t)) / (1024.0 * 1024.0) << " MiB\n";
    cout << " Threads Enabled = " << omp_get_max_threads() << "\n";
    cout << "======================================================================\n";

    auto t0 = high_resolution_clock::now();

    // 2D Array: Rows = Primes, Cols = Sequence terms 0 to R
    vector<vector<uint64_t>> residues(primes.size(), vector<uint64_t>(R + 1, 0));

    // Sequential Prime Evaluation
    for (int i = 0; i < static_cast<int>(primes.size()); ++i) {
        auto t_prime_start = high_resolution_clock::now();
        cout << "\n[PRIME " << i + 1 << "/5] Modulus p = " << primes[i] << " (Primitive Root: " << roots[i] << ")\n";

        compute_9ball_sequence(R, M, primes[i], roots[i], residues[i]);

        auto t_prime_end = high_resolution_clock::now();
        double elapsed_sec = duration_cast<milliseconds>(t_prime_end - t_prime_start).count() / 1000.0;
        cout << "  -> Prime " << i + 1 << " Completed in " << fixed << setprecision(2) << elapsed_sec << " seconds.\n";
    }

    cout << "\nAll transforms complete. Executing 256-bit CRT reconstruction & File I/O...\n";

    string filename = "b055415_NTT.txt";
    ofstream outfile(filename);

    for (uint64_t r = 0; r <= R; ++r) {
        vector<uint64_t> r_slices(primes.size());
        for (size_t p_idx = 0; p_idx < primes.size(); ++p_idx) {
            r_slices[p_idx] = residues[p_idx][r];
        }
        string exact_N9 = solve_crt_garner_256(primes, r_slices);
        outfile << r << " " << exact_N9 << "\n";
    }

    outfile.close();

    auto t1 = high_resolution_clock::now();
    double total_sec = duration_cast<milliseconds>(t1 - t0).count() / 1000.0;

    cout << "----------------------------------------------------------------------\n";
    cout << " Sequence generated successfully -> " << filename << "\n";
    cout << " Total Execution Time: " << fixed << setprecision(3) << total_sec << " seconds.\n";
    cout << "======================================================================\n";

    return 0;
}
