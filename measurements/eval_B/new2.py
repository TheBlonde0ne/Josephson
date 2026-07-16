import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.constants as const

# 1. Load the summary CSV (Must contain the '_Error_uA' columns from the previous step)
df = pd.read_csv('josephson_summary_results.csv')

# 2. Extract B-field from filename
def extract_coil_voltage(filename):
    # e.g., "-0p5V.txt" -> "-0.5"
    base = filename.replace(".txt", "").replace("V", "").replace("p", ".")
    return float(base)

df['Coil_Voltage_V'] = df['File'].apply(extract_coil_voltage)
COIL_CONSTANT = 1.927 / 5.8
df['B_field_mT'] = df['Coil_Voltage_V'] * COIL_CONSTANT

# 3. Filter only superconducting states and calculate absolute average critical current + Errors
df_super = df[df['State'] == 'Superconducting'].copy()

# Average the jump magnitudes to get Ic (converted to mA)
df_super['Ic_mA'] = (np.abs(df_super['Pos_Jump_Current_uA']) + np.abs(df_super['Neg_Jump_Current_uA'])) / 2.0 / 1000.0

# Propagate the polyfit errors in quadrature (converted to mA)
df_super['Ic_Error_mA'] = (0.5 * np.sqrt(
    df_super['Pos_Jump_Current_Error_uA']**2 + 
    df_super['Neg_Jump_Current_Error_uA']**2
)) / 1000.0

# Sort monotonically by magnetic field for clean plotting/fitting
df_super = df_super.sort_values(by='B_field_mT').reset_index(drop=True)

B_mT = df_super['B_field_mT'].values
Ic_mA = df_super['Ic_mA'].values
Ic_err = df_super['Ic_Error_mA'].values

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
        # We pass sigma=Ic_err to weight the fit, ensuring noisy points don't skew the curve
        popt, pcov = curve_fit(fraunhofer_model, B_mT, Ic_mA, p0=p0, bounds=bounds, 
                               sigma=Ic_err, absolute_sigma=True)
        I0_fit, B_offset_fit, delta_B_fit, C_fit = popt
        
        # Calculate theoretical area
        PHI_0 = const.h / (2 * const.e)
        Area_eff = PHI_0 / (delta_B_fit * 1e-3)
        
        # 6. VISUALIZATION
        plt.figure(figsize=(9, 6))
        
        # Plot data with Error Bars (styled to match your reference image)
        plt.errorbar(B_mT, Ic_mA, yerr=Ic_err, fmt='o', 
                     color='red', ecolor='red', capsize=3, elinewidth=1, markersize=4, 
                     label='Measurements', zorder=3)
        
        # Create a smooth line for the theoretical fit
        B_fit_line = np.linspace(np.min(B_mT) - 0.5, np.max(B_mT) + 0.5, 500)
        Ic_fit_line = fraunhofer_model(B_fit_line, *popt)
        
        # Fit line styled as a solid blue line like the image reference
        plt.plot(B_fit_line, Ic_fit_line, color='blue', linewidth=1.5, 
                 label=r'Fit $\mathcal{F}_{a,b}$', zorder=2)

        plt.title('Maximaler Josephson-Strom $I_c$ vs. Magnetische Feldstärke', fontsize=14)
        plt.xlabel('Magnetische Feldstärke $B$ [mT]', fontsize=12)
        plt.ylabel('Maximaler Josephson-Strom $I_c$ [\u03BCA]', fontsize=12) # Kept label unit matching the image
        
        # Use a more subtle layout to match standard physics plots
        plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
        plt.gca().tick_params(direction='in', top=True, right=True)
        
        fit_text = (f"Fit Results:\n"
                    f"$I_0$ = {I0_fit:.4f} mA\n"
                    f"$\\Delta B$ = {delta_B_fit:.4f} mT\n"
                    f"$B_{{offset}}$ = {B_offset_fit:.4f} mT\n"
                    f"$C$ = {C_fit:.4f} mA\n\n"
                    f"$A_{{eff}}$ = {Area_eff:.2e} m$^2$")
        plt.text(0.25, 0.95, fit_text, transform=plt.gca().transAxes, fontsize=11,
                 verticalalignment='top', horizontalalignment='right', 
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
        
        # --- Insert this right after your curve_fit block inside the try-except ---

        # 1. Sample 1000 variations of parameters based on the fit covariance matrix
        # This automatically handles the correlations between your fit parameters
        np.random.seed(42) # For reproducible bands
        psample = np.random.multivariate_normal(popt, pcov, size=1000)

        # 2. Evaluate the model at every point along the smooth B_fit_line for ALL 1000 samples
        # This creates a 2D array of shape (1000, 500)
        sample_curves = np.array([fraunhofer_model(B_fit_line, *p) for p in psample])

        # 3. Calculate the standard deviation (1-sigma range) at each point along the line
        # Axis 0 calculates the variation across the 1000 simulated curves
        sigma_curve = np.std(sample_curves, axis=0)

        # 4. Define the upper and lower bounds for the 1-sigma interval
        Ic_fit_upper = Ic_fit_line + sigma_curve
        Ic_fit_lower = np.maximum(0, Ic_fit_line - sigma_curve) # Avoid physically impossible negative currents

        # 5. Plotting the shaded region
        plt.fill_between(B_fit_line, Ic_fit_lower, Ic_fit_upper, 
                 color='#0047AB', alpha=0.15, label=r'$\sigma_1$ Fit Confidence Band', zorder=1)


        plt.legend(loc='upper right', frameon=False)
        plt.tight_layout()
        plt.savefig("fraunhofer_fit.svg")
        plt.savefig("fraunhofer_fit.png", dpi=300)
        
        # 7. EXPORT DATA (Added Error column to the export file)
        Ic_fit_values = fraunhofer_model(B_mT, *popt)
        np.savetxt("extracted_Ic_vs_B.txt", np.column_stack((B_mT, Ic_mA, Ic_err, Ic_fit_values)), 
                   header="B_mT\tIc_mA\tIc_Err_mA\tFit_mA", fmt="%.6e", comments="")
                   
        np.savetxt("fit_parameters.txt", np.column_stack((I0_fit, B_offset_fit, delta_B_fit, Area_eff)), 
                   header="I0_mA\tB_offset_mT\tdelta_B_mT\tArea_eff_m2", fmt="%.6e", comments="")

        print("\nSuccessfully fit Fraunhofer pattern!")
        print(f"I0 = {I0_fit:.4f} mA, delta_B = {delta_B_fit:.4f} mT")
        
    except Exception as e:
        print(f"Fitting failed: {e}")