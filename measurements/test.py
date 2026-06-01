import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as const

# ==========================================
# 1. PARAMETERS & CONSTANTS
# ==========================================
DATA_DIR = 'Bdependence'     # Folder containing the -0p5V.txt files
V_CRITERION = 4e-5           # 50 uV threshold (applied AFTER offset correction)
SKIP_ROWS = 1                # Skip the text header row

# Helmholtz Coil Calibration (from Section 9.1)
# B/I = (4/5)^(3/2) * mu_0 * (n/R) = 1.927 mT/A
# R_coil approx 5.8 Ohms (from temperatures.txt)
COIL_CONSTANT = 1.927 / 5.8  # Conversion factor: Coil Voltage (V) to B-field (mT)

# Physical constants for theoretical extraction
PHI_0 = const.h / (2 * const.e)
a_length = 10e-6             
D_barrier = 30e-9            

# ==========================================
# 2. DATA EXTRACTION ROUTINE
# ==========================================
def extract_ic_from_vi(voltage_array, current_array, threshold):
    """Finds the current where voltage exceeds the superconducting threshold."""
    # Ensure arrays are sorted by absolute current
    sort_idx = np.argsort((current_array))
    v_sorted = (voltage_array[sort_idx])
    i_sorted = (current_array[sort_idx])
    # sort_idx = np.argsort(np.abs(current_array))
    # v_sorted = np.abs(voltage_array[sort_idx])
    # i_sorted = np.abs(current_array[sort_idx])
    
    # Find where voltage crosses the criterion
    transition_indices = np.where(v_sorted > threshold)[0]

    value = 0

    v_rel = [v_sorted[i] - (i*value) for i in range(len(v_sorted))]

    # plot for debugging
    plt.figure(figsize=(6, 4))
    plt.plot(v_sorted * 1e6, i_sorted * 1e6, marker='o', linestyle='-', color='blue', markersize=1, label='Measured I-V')
    # plt.axhline(threshold * 1e6, color='red', linestyle='--', label='Voltage Criterion')
    plt.ylabel('Current (µA)')
    plt.xlabel('Voltage (µV)')
    plt.title('I-V Curve with Criterion')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    # plt.xscale('log')
    # plt.yscale('log')
    plt.show()
    
    if len(transition_indices) > 0:
        return i_sorted[transition_indices[0]]
    else:
        return np.nan # Return NaN if no transition is found

# Initialize lists to hold our processed data
coil_voltages = []
extracted_Ic = []

# Find all txt files in the directory
file_pattern = os.path.join(DATA_DIR, '*.txt')
files = glob.glob(file_pattern)

if not files:
    print(f"Warning: No files found in directory '{DATA_DIR}'. Check the path.")

for file_path in files:
    # 1. Parse filename (e.g., "-0p5V.txt" -> -0.5)
    filename = os.path.basename(file_path)
    val_str = filename.replace('.txt', '').replace('V', '').replace('p', '.')
    
    try:
        coil_v = float(val_str)
    except ValueError:
        continue # Skip files that don't match the naming convention
        
    # 2. Load V-I data
    try:
        data = np.genfromtxt(file_path, skip_header=SKIP_ROWS, usecols=(0, 1))
        meas_v = data[:, 0]
        meas_i = data[:, 1]
        
        # --- CRITICAL FIX: OFFSET REMOVAL ---
        # Find the baseline voltage around zero current (e.g., between -5uA and +5uA)
        zero_current_mask = np.abs(meas_i) < 5e-6
        if np.any(zero_current_mask):
            v_offset = np.mean(meas_v[zero_current_mask])
        else:
            v_offset = 0.0
            
        meas_v_corrected = meas_v - v_offset
        # ------------------------------------
        
        # 3. Extract critical current using corrected voltage
        ic = extract_ic_from_vi(meas_v_corrected, meas_i, V_CRITERION)
        
        if not np.isnan(ic):
            coil_voltages.append(coil_v)
            extracted_Ic.append(ic)
        
    except Exception as e:
        print(f"Error reading {filename}: {e}")

# Convert to numpy arrays and sort by coil voltage
coil_voltages = np.array(coil_voltages)
extracted_Ic = np.array(extracted_Ic)

sort_idx = np.argsort(coil_voltages)
coil_voltages = coil_voltages[sort_idx]
extracted_Ic = extracted_Ic[sort_idx]

# Convert parameters for fitting
B_mT = coil_voltages * COIL_CONSTANT
Ic_mA = extracted_Ic * 1000  # Convert Amps to milliAmps

# ==========================================
# 3. THEORETICAL MODEL & FITTING
# ==========================================
def fraunhofer_model(B, I0, B_offset, delta_B):
    x = (B - B_offset) / delta_B
    return I0 * np.abs(np.sinc(x))

if len(B_mT) > 0:
    I0_guess = np.max(Ic_mA)
    B_offset_guess = B_mT[np.argmax(Ic_mA)]
    delta_B_guess = (np.max(B_mT) - np.min(B_mT)) / 4.0 

    p0 = [I0_guess, B_offset_guess, delta_B_guess]
    bounds = ([0, -np.inf, 0.001], [np.inf, np.inf, np.inf])

    try:
        popt, pcov = curve_fit(fraunhofer_model, B_mT, Ic_mA, p0=p0, bounds=bounds)
        I0_fit, B_offset_fit, delta_B_fit = popt
        
        # ==========================================
        # 4. VISUALIZATION
        # ==========================================
        plt.figure(figsize=(9, 6))
        plt.scatter(B_mT, Ic_mA, color='red', s=40, label='Extracted $I_c$', zorder=3)
        
        B_fit_line = np.linspace(np.min(B_mT), np.max(B_mT), 500)
        Ic_fit_line = fraunhofer_model(B_fit_line, *popt)
        plt.plot(B_fit_line, Ic_fit_line, color='blue', linewidth=2, label='Fraunhofer Fit')

        plt.title('Magnetic Field Dependence of Critical Current', fontsize=14)
        plt.xlabel('Magnetic Field $B$ (mT)', fontsize=12)
        plt.ylabel('Critical Current $I_c$ (mA)', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        fit_text = (f"Fit Results:\n"
                    f"$I_0$ = {I0_fit:.3f} mA\n"
                    f"$\Delta B$ = {delta_B_fit:.3f} mT\n"
                    f"$B_{{offset}}$ = {B_offset_fit:.3f} mT")
        plt.text(0.95, 0.95, fit_text, transform=plt.gca().transAxes, fontsize=11,
                 verticalalignment='top', horizontalalignment='right', 
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
                 
        plt.legend()
        plt.tight_layout()
        plt.savefig("critical_current_fit.svg")

        # Save extracted data and fit parameters

        Ic_fit_values = fraunhofer_model(B_mT, *popt)

        np.savetxt("extracted_Ic_vs_B.txt", np.column_stack((B_mT, Ic_mA, Ic_fit_values)), header="B_mT\tIc_mA\tFit_mA", fmt="%.6e")
        np.savetxt("fit_parameters.txt", np.array(popt), header="I0_mA\tB_offset_mT\tdelta_B_mT", fmt="%.6e")


    except Exception as e:
        print(f"Fitting failed: {e}")