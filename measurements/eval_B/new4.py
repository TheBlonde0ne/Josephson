import os
import numpy as np
import scipy.constants as const

# ==========================================
# 1. PARAMETER-DEKLARATION & DATEIEINGABE
# ==========================================

# Standardwerte (aus Ihrem Canvas-Fit als Fallback)
I0_fit = 5.909298e-02          # mA
I0_err = 3.525113e-03          # mA
delta_B_fit = 1.090710e+00     # mT
delta_B_err = 1.909656e-02     # mT
Area_eff = 1.895861e-12        # m^2
Area_eff_err = 3.319347e-14    # m^2

# Geometrische Abmessungen der Probe (Anpassung je nach Versuch)
w_width = 10.0e-6              # Breite w = 10 um (in Metern)
L_length = 10.0e-6             # Länge L = 10 um (in Metern)
D_barrier = 30.0e-9            # Barrierendicke D = 30 nm (in Metern)

# Versuchen, die Parameter direkt aus der zuvor erstellten Datei zu lesen
param_file = "fit_parameters.txt"
if os.path.exists(param_file):
    try:
        # Einlesen der Kopfzeile und der Werte
        data = np.loadtxt(param_file, skiprows=1)
        if data.ndim == 1 and len(data) >= 10:
            I0_fit, I0_err = data[0], data[1]
            # delta_B an Index 4 und 5
            delta_B_fit, delta_B_err = data[4], data[5]
            # Area_eff an Index 8 und 9
            Area_eff, Area_eff_err = data[8], data[9]
            print(f"-> Parameter erfolgreich aus '{param_file}' geladen.")
    except Exception as e:
        print(f"-> Fehler beim Lesen von '{param_file}': {e}. Nutze Standard-Fallback-Werte.")
else:
    print(f"-> '{param_file}' nicht gefunden. Nutze voreingestellte Standardwerte.")

# Konvertierung in SI-Einheiten
I0_A = I0_fit * 1e-3
I0_err_A = I0_err * 1e-3

# ==========================================
# 2. PHYSIKALISCHE FORMELN & FEHLERFORTPFLANZUNG
# ==========================================

# Konstanten
Phi_0 = const.h / (2 * const.e)  # Magnetisches Flussquant (~ 2.0678e-15 Wb)
mu_0 = const.mu_0                # Vakuumpermeabilität (~ 1.2566e-6 H/m)

# --- A) London'sche Eindringtiefe (lambda_L) ---
# Unter der Annahme symmetrischer Elektroden gilt: A_eff = w * (D + 2*lambda_L)
# Auflösen nach lambda_L:
lambda_L = 0.5 * ((Area_eff / w_width) - D_barrier)

# Fehlerfortpflanzung für lambda_L (nur abhängig vom Fehler der effektiven Fläche Area_eff):
# d(lambda_L)/d(Area_eff) = 1 / (2 * w)
lambda_L_err = Area_eff_err / (2.0 * w_width)


# --- B) Josephson-Eindringtiefe (lambda_J) ---
# Die magnetische Dicke d_mag ist definiert als d_mag = D + 2*lambda_L = A_eff / w
d_mag = D_barrier + 2.0 * lambda_L

# Kritische Stromdichte j_c = I_0 / (w * L)
j_c = I0_A / (w_width * L_length)
# Fehler der kritischen Stromdichte:
j_c_err = I0_err_A / (w_width * L_length)

# Berechnung von lambda_J:
# lambda_J = sqrt( Phi_0 / (2 * pi * mu_0 * d_mag * j_c) )
# Einsetzen von d_mag und j_c zeigt, dass lambda_J = w * sqrt( Phi_0 * L / (2 * pi * mu_0 * Area_eff * I_0) )
denominator = 2.0 * np.pi * mu_0 * d_mag * j_c
lambda_J = np.sqrt(Phi_0 / denominator)

# Fehlerfortpflanzung für lambda_J (abhängig von den unabhängigen Variablen Area_eff und I_0):
# Der relative Fehler lässt sich elegant vereinfachen zu:
# rel_err_lambda_J = 0.5 * sqrt( (Area_err/Area_eff)^2 + (I0_err/I0)^2 )
rel_err_lambda_J = 0.5 * np.sqrt((Area_eff_err / Area_eff)**2 + (I0_err_A / I0_A)**2)
lambda_J_err = lambda_J * rel_err_lambda_J


# ==========================================
# 3. AUSGABE & EXPORT DER ERGEBNISSE
# ==========================================

print("\n" + "="*55)
print(" PHYSIKALISCHE ERGEBNISSE & EINDRINGTIEFEN")
print("="*55)
print(f"Effektive Barrierefläche A_eff:  ({Area_eff*1e12:.4f} \u00b1 {Area_eff_err*1e12:.4f}) \u03bcm\u00b2")
print(f"Maximaler Josephson-Strom I_0:   ({I0_fit*1e3:.2f} \u00b1 {I0_err*1e3:.2f}) \u03bcA")
print(f"Kritische Stromdichte j_c:       ({j_c/1e2:.2f} \u00b1 {j_c_err/1e2:.2f}) A/cm\u00b2")
print("-"*55)
print(f"London'sche Eindringtiefe \u03bb_L:   ({lambda_L*1e9:.2f} \u00b1 {lambda_L_err*1e9:.2f}) nm")
print(f"Josephson-Eindringtiefe \u03bb_J:    ({lambda_J*1e6:.2f} \u00b1 {lambda_J_err*1e6:.2f}) \u03bcm")
print("="*55)

# Verhältnis Junction-Länge L zu lambda_J bewerten ("kurzer" vs "langer" Kontakt)
verhaeltnis = L_length / lambda_J
print(f"Verhältnis L / \u03bb_J:               {verhaeltnis:.4f}")
if verhaeltnis < 1.0:
    print("-> Da L < \u03bb_J, handelt es sich um einen 'kurzen' Josephson-Kontakt.")
    print("   Die Annahme eines homogenen Magnetfeldes im Kontakt ist gerechtfertigt.")
else:
    print("-> Da L >= \u03bb_J, handelt es sich um einen 'langen' Josephson-Kontakt.")
    print("   Es treten Selbstabschirmungseffekte und Josephson-Wirbel auf.")
print("="*55)

# Export in Datei für das Protokoll
export_file = "eindringtiefen_ergebnisse.txt"
with open(export_file, "w", encoding="utf-8") as f:
    f.write("Ergebnisse der physikalischen Eindringtiefen-Analyse\n")
    f.write("===================================================\n\n")
    f.write(f"London'sche Eindringtiefe lambda_L:   ({lambda_L*1e9:.4f} +/- {lambda_L_err*1e9:.4f}) nm\n")
    f.write(f"Josephson-Eindringtiefe lambda_J:    ({lambda_J*1e6:.4f} +/- {lambda_J_err*1e6:.4f}) um\n")
    f.write(f"Kritische Stromdichte j_c:           ({j_c/1e2:.4f} +/- {j_c_err/1e2:.4f}) A/cm^2\n")
    f.write(f"Verhaeltnis L / lambda_J:            {verhaeltnis:.4f}\n")

print(f"\nErgebnisse wurden erfolgreich in '{export_file}' exportiert.")