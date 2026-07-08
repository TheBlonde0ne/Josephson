import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def get_fit_errors(cov):
    """
    Helper to extract slope/intercept errors. 
    If cov is a float (SEM), it assumes slope error is 0 (constant line).
    """
    if isinstance(cov, (float, int, np.float64)):
        return 0.0, float(cov)
    else:
        return np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1])

def find_intersection_with_error(I, fit1, cov1, fit2, cov2):
    """
    Finds intersection and propagates standard errors from polyfit.
    """
    a1, b1 = fit1[0], fit1[1]
    a2, b2 = fit2[0], fit2[1]
    
    # Extract standard errors using the helper
    sig_a1, sig_b1 = get_fit_errors(cov1)
    sig_a2, sig_b2 = get_fit_errors(cov2)
        
    # Calculate intersection
    I_int = (b2 - b1) / (a1 - a2)
    V_int = a1 * I_int + b1
    
    # Partial derivatives for error propagation
    dI_da1 = (b1 - b2) / (a1 - a2)**2
    dI_da2 = (b2 - b1) / (a1 - a2)**2
    dI_db1 = -1 / (a1 - a2)
    dI_db2 = 1 / (a1 - a2)
    
    # Propagate variance
    I_var = (dI_da1 * sig_a1)**2 + (dI_db1 * sig_b1)**2 + (dI_da2 * sig_a2)**2 + (dI_db2 * sig_b2)**2
    sig_I_int = np.sqrt(I_var)
    # Calculate the mean spacing of your current array as a resolution floor
    I_resolution_floor = np.mean(np.diff(I)) 

    # If the propagated error explodes due to steep slopes, cap it reasonably
    # A jump error shouldn't realistically exceed 2-3 measurement steps
    max_allowable_error = 3 * I_resolution_floor 

    # if sig_I_int > max_allowable_error:
    #     sig_I_int = I_resolution_floor  # Fall back to your instrument's step size resolution
    return I_int, V_int, sig_I_int

def check_resistive_state(I, V, r_squared_threshold=0.9993):
    """
    Checks if the entire dataset is a single straight line (no superconducting gap).
    Returns True if pure Ohmic resistance, False if Josephson jumps are present.
    """
    fit = np.polyfit(I, V, 1)
    V_pred = np.polyval(fit, I)
    
    ss_res = np.sum((V - V_pred)**2)
    ss_tot = np.sum((V - np.mean(V))**2)
    r2 = 1 - (ss_res / ss_tot)
    
    return (r2) > r_squared_threshold, fit[0]

def evaluate_josephson_file(file_path):
    df = pd.read_csv(file_path, sep="\t")
    v_col, i_col = df.columns[0], df.columns[1]
    
    df = df.sort_values(by=i_col).reset_index(drop=True)
    V = df[v_col].values
    I = df[i_col].values
    
    # ------------------------------------------------------------------
    # Pre-Check: Is the junction fully resistive? (No Superconductivity)
    # ------------------------------------------------------------------
    is_resistive, normal_R = check_resistive_state(I, V)
    filename = os.path.basename(file_path)
    
    if is_resistive:
        # Plot the single straight line
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(I * 1e6, V * 1e3, 'k.', markersize=2, alpha=0.5, label='Data (Fully Resistive)')
        i_line = np.linspace(I.min(), I.max(), 100)
        ax.plot(i_line * 1e6, (normal_R * i_line + np.polyfit(I, V, 1)[1]) * 1e3, 'r-', linewidth=2, label=f'Ohmic Fit (R={normal_R:.2f}\u03A9)')
        ax.set_title(f'No Superconducting State Detected: {filename}')
        ax.set_xlabel('Current (\u03bcA)')
        ax.set_ylabel('Voltage (mV)')
        ax.grid(True, linestyle=':', alpha=0.6)
        ax.legend()
        
        os.makedirs("plots", exist_ok=True)
        plt.savefig(f"plots/fit_{filename.replace('.txt', '.png')}", dpi=200)
        plt.close()
        
        return {
            "File": filename,
            "State": "Resistive",
            "Pos_Jump_Current_uA": 0,
            "Pos_Jump_Current_Error_uA": 0,
            "Neg_Jump_Current_uA": 0,
            "Neg_Jump_Current_Error_uA": 0,
            "Pos_Jump_Voltage_mV": 0,
            "Neg_Jump_Voltage_mV": 0,
            "Plateau_Offset_mV": 0,
            "Normal_Resistance_Ohms": normal_R
        }

    # ------------------------------------------------------------------
    # Standard Analysis for Superconducting State
    # ------------------------------------------------------------------
    dV = np.gradient(V)
    dI = np.gradient(I)
    
    valid = np.abs(dI) > 1e-12
    R = np.zeros_like(V)
    R[valid] = np.abs(dV[valid] / dI[valid])
    
    neg_jump_center = np.argmax(R * (I < 0))
    pos_jump_center = np.argmax(R * (I >= 0))
    
    R_branch = np.median(R[0:min(50, len(R)//10)])
    threshold = max(2.0 * R_branch, 10.0)
    
    j2_start, j2_end = neg_jump_center, neg_jump_center
    while j2_start > 0 and R[j2_start] > threshold: j2_start -= 1
    while j2_end < len(R)-1 and R[j2_end] > threshold: j2_end += 1
        
    j1_start, j1_end = pos_jump_center, pos_jump_center
    while j1_start > 0 and R[j1_start] > threshold: j1_start -= 1
    while j1_end < len(R)-1 and R[j1_end] > threshold: j1_end += 1
        
    buf = 2
    idx_b2   = range(0, max(0, j2_start - buf))
    idx_j2   = range(j2_start, j2_end + 1)
    idx_plat = range(j2_end + buf, max(j2_end + buf + 1, j1_start - buf))
    idx_j1   = range(j1_start, j1_end + 1)
    idx_b1   = range(j1_end + buf, len(I))
    
    # Generate Fits and Covariances
    fit_b2, cov_b2     = np.polyfit(I[idx_b2], V[idx_b2], 1, cov=True)
    fit_j2, cov_j2     = np.polyfit(I[idx_j2], V[idx_j2], 1, cov=True)
    fit_plat, cov_plat = np.polyfit(I[idx_plat], V[idx_plat], 1, cov=True)
    fit_j1, cov_j1     = np.polyfit(I[idx_j1], V[idx_j1], 1, cov=True)
    fit_b1, cov_b1     = np.polyfit(I[idx_b1], V[idx_b1], 1, cov=True)

    mean_plat = np.mean(V[idx_plat])
    cov_plat = np.std(V[idx_plat], ddof=1) / np.sqrt(len(idx_plat)) if len(idx_plat) > 1 else 1e-6
    fit_plat = np.array([0.0, mean_plat])
    
    # Calculate Intersections with Error Propagation
    I_neg_top, V_neg_top, sig_neg_top = find_intersection_with_error(I, fit_b2, cov_b2, fit_j2, cov_j2)
    I_neg_bot, V_neg_bot, sig_neg_bot = find_intersection_with_error(I, fit_j2, cov_j2, fit_plat, cov_plat)
    I_pos_bot, V_pos_bot, sig_pos_bot = find_intersection_with_error(I, fit_plat, cov_plat, fit_j1, cov_j1)
    I_pos_top, V_pos_top, sig_pos_top = find_intersection_with_error(I, fit_j1, cov_j1, fit_b1, cov_b1)
    
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(I * 1e6, V * 1e3, 'k.', markersize=2, alpha=0.5, label='Data')
    
    i_ranges = [
        (I[idx_b2][0], I[idx_b2][-1]), (I[idx_j2][0], I[idx_j2][-1]), 
        (I[idx_plat][0], I[idx_plat][-1]), (I[idx_j1][0], I[idx_j1][-1]), 
        (I[idx_b1][0], I[idx_b1][-1])
    ]
    fits = [fit_b2, fit_j2, fit_plat, fit_j1, fit_b1]
    colors = ['c-', 'm--', 'b-', 'g--', 'r-']
    
    for (i_min, i_max), fit, col in zip(i_ranges, fits, colors):
        i_line = np.linspace(i_min, i_max, 50)
        ax.plot(i_line * 1e6, (fit[0] * i_line + fit[1]) * 1e3, col, linewidth=2)

    ints = [(I_neg_top, V_neg_top, 'mo'), (I_neg_bot, V_neg_bot, 'bo'), 
            (I_pos_bot, V_pos_bot, 'go'), (I_pos_top, V_pos_top, 'ro')]
    for int_I, int_V, m in ints:
        ax.plot(int_I * 1e6, int_V * 1e3, m, markersize=6)
        
    ax.set_title(f'Superconducting State Identified: {filename}')
    ax.set_xlabel('Current (\u03bcA)')
    ax.set_ylabel('Voltage (mV)')
    ax.grid(True, linestyle=':', alpha=0.6)
    
    os.makedirs("plots", exist_ok=True)
    plt.savefig(f"plots/fit_{filename.replace('.txt', '.png')}", dpi=200)
    plt.close()
    
    return {
        "File": filename,
        "State": "Superconducting",
        "Pos_Jump_Current_uA": I_pos_bot * 1e6,
        "Pos_Jump_Current_Error_uA": sig_pos_bot * 1e6,
        "Neg_Jump_Current_uA": I_neg_bot * 1e6,
        "Neg_Jump_Current_Error_uA": sig_neg_bot * 1e6,
        "Pos_Jump_Voltage_mV": V_pos_bot * 1e3,
        "Neg_Jump_Voltage_mV": V_neg_bot * 1e3,
        "Plateau_Offset_mV": fit_plat[1] * 1e3,
        "Normal_Resistance_Ohms": (fit_b1[0] + fit_b2[0]) / 2
    }

def process_all_files(directory="C:\\Users\\walln\\Documents\\Uni\\Josephson\\measurements\\Bdependence"):
    data_files = glob.glob(os.path.join(directory, "*.txt"))
    results = []
    
    for f in data_files:
        try:
            file_results = evaluate_josephson_file(f)
            results.append(file_results)
            print(f"[{file_results['State']}] Processed {os.path.basename(f)}")
        except Exception as e:
            print(f"[ERROR] Failed processing {f}: {e}")
            
    if results:
        pd.DataFrame(results).to_csv("josephson_summary_results.csv", index=False)
        print("\nSummary saved to 'josephson_summary_results.csv'.")

# Run the batch processor
if __name__ == "__main__":
    process_all_files("C:\\Users\\walln\\Documents\\Uni\\Josephson\\measurements\\Bdependence")