import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as const

# 1. Load the summary CSV (Assuming it's named 'josephson_summary_results.csv')
df = pd.read_csv('josephson_summary_results.csv')

# 2. Extract B-field from filename
def extract_coil_voltage(filename):
    # e.g., "-0p5V.txt" -> "-0.5"
    base = filename.replace(".txt", "").replace("V", "").replace("p", ".")
    return float(base)

df['Coil_Voltage_V'] = df['File'].apply(extract_coil_voltage)
COIL_CONSTANT = 1.927 / 5.8
df['B_field_mT'] = df['Coil_Voltage_V'] * COIL_CONSTANT

# 3. Filter only superconducting states and calculate absolute average critical current
df_super = df[df['State'] == 'Superconducting'].copy()
df_super['Ic_mA'] = (np.abs(df_super['Pos_Jump_Current_uA']) + np.abs(df_super['Neg_Jump_Current_uA'])) / 2.0 / 1000.0

# Sort monotonically by magnetic field for clean plotting/fitting
df_super = df_super.sort_values(by='B_field_mT').reset_index(drop=True)

B_mT = df_super['B_field_mT'].values
Ic_mA = df_super['Ic_mA'].values

# 4. Define Fraunhofer Physics Model
def fraunhofer_model(B, I0, B_offset, delta_B, C):
    # np.sinc(x) is sin(pi*x)/(pi*x)
    x = (B - B_offset) / delta_B
    return I0 * np.abs(np.sinc(x)) + C

if len(B_mT) > 0:
    # 5. Smart Initial Guesses
    I0_guess = np.max(Ic_mA)
    B_offset_guess = B_mT[np.argmax(Ic_mA)]
    delta_B_guess = (np.max(B_mT) - np.min(B_mT)) / 4.0 
    C_guess = np.min(Ic_mA)

    p0 = [I0_guess, B_offset_guess, delta_B_guess, C_guess]
    bounds = ([0, -np.inf, 0.001, -np.inf], [np.inf, np.inf, np.inf, np.inf])

    try:
        # Execute non-linear least squares fit
        popt, pcov = curve_fit(fraunhofer_model, B_mT, Ic_mA, p0=p0, bounds=bounds)
        I0_fit, B_offset_fit, delta_B_fit, C_fit = popt
        
        # Calculate theoretical area
        PHI_0 = const.h / (2 * const.e)
        Area_eff = PHI_0 / (delta_B_fit * 1e-3)
        
        # 6. VISUALIZATION
        plt.figure(figsize=(9, 6))
        plt.scatter(B_mT, Ic_mA, edgecolors='#F54927', color="white", s=45, label='Extracted $I_c$', zorder=3)
        
        # Create a smooth line for the theoretical fit
        B_fit_line = np.linspace(np.min(B_mT) - 0.5, np.max(B_mT) + 0.5, 500)
        Ic_fit_line = fraunhofer_model(B_fit_line, *popt)
        plt.plot(B_fit_line, Ic_fit_line, color='#0047AB', linewidth=2, linestyle="--", label='Fraunhofer Fit')

        plt.title('Josephson Critical Current vs. Magnetic Field', fontsize=14, fontweight='bold')
        plt.xlabel('Magnetic Field $B$ (mT)', fontsize=12)
        plt.ylabel('Critical Current $I_c$ (mA)', fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        
        fit_text = (f"Fit Results:\n"
                    f"$I_0$ = {I0_fit:.4f} mA\n"
                    f"$\\Delta B$ = {delta_B_fit:.4f} mT\n"
                    f"$B_{{offset}}$ = {B_offset_fit:.4f} mT\n"
                    f"$C$ = {C_fit:.4f} mA\n\n"
                    f"$A_{{eff}}$ = {Area_eff:.2e} m$^2$")
        plt.text(0.95, 0.95, fit_text, transform=plt.gca().transAxes, fontsize=11,
                 verticalalignment='top', horizontalalignment='right', 
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
                 
        plt.legend(loc='upper left')
        plt.tight_layout()
        plt.savefig("fraunhofer_fit.svg")
        plt.savefig("fraunhofer_fit.png", dpi=300)
        
        # 7. EXPORT DATA
        Ic_fit_values = fraunhofer_model(B_mT, *popt)
        np.savetxt("extracted_Ic_vs_B.txt", np.column_stack((B_mT, Ic_mA, Ic_fit_values)), 
                   header="B_mT\tIc_mA\tFit_mA", fmt="%.6e", comments="")
                   
        np.savetxt("fit_parameters.txt", np.column_stack((I0_fit, B_offset_fit, delta_B_fit, Area_eff)), 
                   header="I0_mA\tB_offset_mT\tdelta_B_mT\tArea_eff_m2", fmt="%.6e", comments="")

        print("\nSuccessfully fit Fraunhofer pattern!")
        print(f"I0 = {I0_fit:.4f} mA, delta_B = {delta_B_fit:.4f} mT")
        
    except Exception as e:
        print(f"Fitting failed: {e}")