import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

path = "C:/Users/walln/Documents/Uni/Josephson/measurements/Bdependence/-0V.txt"

df = pd.read_csv(path, sep="\t")
v_col, i_col = df.columns[0], df.columns[1]
df = df.sort_values(by=i_col).reset_index(drop=True)
V = df[v_col].values
I = df[i_col].values

plt.figure(figsize=(8, 5))
plt.grid(which="both", alpha=0.5, linestyle=":")
plt.scatter(I * 1e6, V * 1e3, s=10, edgecolor="#002D6D", label='Data Points')
plt.legend()
plt.xlabel('Current (\u03bcA)')
plt.ylabel('Voltage (mV)')
plt.savefig("simple_curve.png", dpi=200)
plt.close()