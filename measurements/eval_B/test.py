import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def find_intersection(fit1, fit2):
    """Calculates the analytical intersection (I, V) of two linear fits."""
    a1, b1 = fit1[0], fit1[1]
    a2, b2 = fit2[0], fit2[1]
    I_int = (b2 - b1) / (a1 - a2)
    V_int = a1 * I_int + b1
    return I_int, V_int

def evaluate_josephson_file(file_path):
    """Analyzes a single I-V data file and dynamically identifies switching regimes."""
    # Load data (Adjust sep="\t" or sep="," based on your actual files)
    df = pd.read_csv(file_path, sep="\t")
    
    # Assume Column 0 is Voltage and Column 1 is Current
    v_col, i_col = df.columns[0], df.columns[1]
    
    # 1. Sort by Current (ascending) ensures consistent algorithmic logic
    df = df.sort_values(by=i_col).reset_index(drop=True)
    V = df[v_col].values
    I = df[i_col].values
    
    # 2. Compute absolute derivative R = |dV / dI|
    dV = np.gradient(V)
    dI = np.gradient(I)
    
    # Safely handle potential division by zero
    valid = np.abs(dI) > 1e-12
    R = np.zeros_like(V)
    R[valid] = np.abs(dV[valid] / dI[valid])
    
    # 3. Locate Jump Centers (Max Resistance)
    # The negative jump happens when I < 0, positive when I > 0
    neg_jump_center = np.argmax(R * (I < 0))
    pos_jump_center = np.argmax(R * (I >= 0))
    
    # 4. Dynamically expand outwards to find the full jump boundary
    # Estimate the baseline branch resistance from the first 50 points
    R_branch = np.median(R[0:min(50, len(R)//10)])
    threshold = max(2.0 * R_branch, 10.0) # Jump threshold condition
    
    j2_start, j2_end = neg_jump_center, neg_jump_center
    while j2_start > 0 and R[j2_start] > threshold: j2_start -= 1
    while j2_end < len(R)-1 and R[j2_end] > threshold: j2_end += 1
        
    j1_start, j1_end = pos_jump_center, pos_jump_center
    while j1_start > 0 and R[j1_start] > threshold: j1_start -= 1
    while j1_end < len(R)-1 and R[j1_end] > threshold: j1_end += 1
        
    # 5. Define Safe Index Segments (with a buffer of 2 to avoid rounding corners)
    buf = 2
    idx_b2   = range(0, max(0, j2_start - buf))
    idx_j2   = range(j2_start, j2_end + 1)
    idx_plat = range(j2_end + buf, max(j2_end + buf + 1, j1_start - buf))
    idx_j1   = range(j1_start, j1_end + 1)
    idx_b1   = range(j1_end + buf, len(I))
    
    # 6. Apply Polynomial Fits
    fit_b2   = np.polyfit(I[idx_b2], V[idx_b2], 1)
    fit_j2   = np.polyfit(I[idx_j2], V[idx_j2], 1)
    fit_plat = np.polyfit(I[idx_plat], V[idx_plat], 1)
    fit_j1   = np.polyfit(I[idx_j1], V[idx_j1], 1)
    fit_b1   = np.polyfit(I[idx_b1], V[idx_b1], 1)
    
    # 7. Extract Analytical Intersections
    I_neg_top, V_neg_top = find_intersection(fit_b2, fit_j2)
    I_neg_bot, V_neg_bot = find_intersection(fit_j2, fit_plat)
    I_pos_bot, V_pos_bot = find_intersection(fit_plat, fit_j1)
    I_pos_top, V_pos_top = find_intersection(fit_j1, fit_b1)
    
    # --- VISUALIZATION VERIFICATION ---
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(I * 1e6, V * 1e3, 'k.', markersize=2, alpha=0.5, label='Data')
    
    # Plotting Lines
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

    # Plot Intersections
    ints = [(I_neg_top, V_neg_top, 'mo'), (I_neg_bot, V_neg_bot, 'bo'), 
            (I_pos_bot, V_pos_bot, 'go'), (I_pos_top, V_pos_top, 'ro')]
    for int_I, int_V, m in ints:
        ax.plot(int_I * 1e6, int_V * 1e3, m, markersize=6)
        
    filename = os.path.basename(file_path)
    ax.set_title(f'Dynamic Fit Verification: {filename}')
    ax.set_xlabel('Current (\u03bcA)')
    ax.set_ylabel('Voltage (mV)')
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # Save the plot in a 'plots' folder
    os.makedirs("plots", exist_ok=True)
    plt.savefig(f"plots/fit_{filename.replace('.txt', '.png')}", dpi=200)
    plt.close()
    
    # Return extracted physics parameters
    return {
        "File": filename,
        "Pos_Jump_Current_uA": I_pos_bot * 1e6,
        "Pos_Jump_Voltage_mV": V_pos_bot * 1e3,
        "Neg_Jump_Current_uA": I_neg_bot * 1e6,
        "Neg_Jump_Voltage_mV": V_neg_bot * 1e3,
        "Plateau_Offset_mV": fit_plat[1] * 1e3,
        "Normal_Resistance_Ohms": (fit_b1[0] + fit_b2[0]) / 2 # Average of the two resistive branches
    }

# --- BATCH EXECUTION ---
def process_all_files(directory="C:\\Users\\walln\\Documents\\Uni\\Josephson\\measurements\\Bdependence"):
    print("Starting automated evaluation...")
    data_files = glob.glob(os.path.join(directory, "*.txt")) # Change extension if needed
    
    results = []
    for f in data_files:
        try:
            print(f"Processing {os.path.basename(f)}...")
            file_results = evaluate_josephson_file(f)
            results.append(file_results)
        except Exception as e:
            print(f"Error processing {f}: {e}")
            
    if results:
        # Save all results to a single summary CSV
        summary_df = pd.DataFrame(results)
        summary_df.to_csv("josephson_summary_results.csv", index=False)
        print("\nAll files processed successfully!")
        print("Detailed metrics saved to 'josephson_summary_results.csv'.")
        print("Verification plots saved in the '/plots' directory.")
    else:
        print("No files found or all failed.")

# Run the batch processor
if __name__ == "__main__":
    process_all_files("C:\\Users\\walln\\Documents\\Uni\\Josephson\\measurements\\Bdependence") # <-- Provide your actual folder path here