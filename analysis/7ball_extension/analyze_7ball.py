import numpy as np

def load_b_file(filepath):
    radii = []
    counts = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                radii.append(int(parts[0]))
                counts.append(int(parts[1]))
    return np.array(radii, dtype=np.float64), np.array(counts, dtype=np.float64)

def comprehensive_lattice_analysis(filepath='b055413.txt'):
    print("==================================================")
    print("   ADVANCED 7-BALL LATTICE ERROR SUITE (R <= 5000)  ")
    print("==================================================")
    
    # 1. Load Data
    R, N_R = load_b_file(filepath)
    mask = R > 0
    R = R[mask]
    N_R = N_R[mask]
    
    # 2. Continuous Volume & Error Term
    # V_7(R) = (16 * pi^3 / 105) * R^7
    V_7 = (16.0 * (np.pi**3) / 105.0) * (R**7)
    E_R = N_R - V_7
    
    # 3. Shell Counts (First Differences: points added per integer shell step)
    shell_counts = np.diff(N_R, prepend=N_R[0])
    
    print(f"\n[+] Dataset Loaded: {len(R)} radial steps processed.")
    print(f"[+] Max Radius: {int(R[-1])} | Max N(R): {int(N_R[-1])}")
    
    # 4. Rolling / Bracketed Exponent Analysis (Local Scaling of |E(R)|)
    print("\n--- Local Bounding Exponent (alpha) Across Brackets ---")
    brackets = [(100, 1000), (1000, 2000), (2000, 3000), (3000, 4000), (4000, 5000)]
    for low, high in brackets:
        b_mask = (R >= low) & (R <= high)
        if np.any(b_mask):
            r_sub = R[b_mask]
            e_sub = np.abs(E_R[b_mask])
            log_r = np.log(r_sub)
            log_e = np.log(np.maximum(e_sub, 1e-12))
            alpha_local, _ = np.polyfit(log_r, log_e, 1)
            print(f"  Radius [{low:4d} - {high:4d}] -> Local Exponent alpha: {alpha_local:.4f}")

    # 5. Frequency Domain Analysis (FFT on Error Term for High R)
    # Focus on the high-end extension where data is pristine and unexamined
    high_ext_mask = R >= 1000
    r_ext = R[high_ext_mask]
    e_ext = E_R[high_ext_mask]
    
        # Detrend error to look purely at fluctuations
    poly_detrend = np.polyfit(r_ext, e_ext, 2)
    trend = np.polyval(poly_detrend, r_ext)
    detrended_error = e_ext - trend
    
    # Compute FFT
    fft_vals = np.fft.rfft(detrended_error)
    fft_freqs = np.fft.rfftfreq(len(detrended_error))
    power_spectrum = np.abs(fft_vals)**2
    
    # Find dominant frequency peaks
    peak_idx = np.argsort(power_spectrum)[-3:][::-1] # Top 3 peaks
    print("\n--- Dominant Frequency Modes in High-End Extension (R >= 1000) ---")
    for idx in peak_idx:
        freq = fft_freqs[idx]
        period = 1.0 / freq if freq > 0 else float('inf')
        print(f"  Peak Frequency: {freq:.5f} (Corresponding Period Length: ~{period:.2f} radial units, Power: {power_spectrum[idx]:.2e})")

    # 6. Statistical Profiling of High-End Residuals
    mean_err = np.mean(e_ext)
    std_err = np.std(e_ext)
    print("\n--- High-End Error Distribution Statistics (R >= 1000) ---")
    print(f"  Mean Error: {mean_err:,.2f}")
    print(f"  Standard Deviation: {std_err:,.2f}")
    print(f"  Max Absolute Error: {np.max(np.abs(e_ext)):,.2f}")

    print("\n==================================================")
    print("   ANALYSIS COMPLETE. READY FOR THEORETICAL PARSING. ")
    print("==================================================")

if __name__ == "__main__":
    comprehensive_lattice_analysis()