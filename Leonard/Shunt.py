from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATA_FILE = Path(__file__).with_name("Tdependence") / "4p191K.txt"
OUTPUT_PLOT = Path(__file__).with_name("Plots") / "Shunt_fits_4p191K.pdf"


def load_data(path: Path) -> tuple[np.ndarray, np.ndarray]:
	data = np.loadtxt(path)
	voltage = data[:, 0]
	current = data[:, 1]
	return current, voltage


def fit_linear(current: np.ndarray, voltage: np.ndarray) -> tuple[float, float]:
	coefficients, covariance = np.polyfit(current, voltage, 1, cov=True)
	slope, intercept = coefficients
	slope_error = float(np.sqrt(covariance[0, 0]))
	intercept_error = float(np.sqrt(covariance[1, 1]))
	return float(slope), float(intercept), slope_error, intercept_error


def linear_curve(current: np.ndarray, slope: float, intercept: float) -> np.ndarray:
	return slope * current + intercept


def fit_region(
	current: np.ndarray,
	voltage: np.ndarray,
	*,
	lower_current: float | None = None,
	upper_current: float | None = None,
	lower_voltage: float | None = None,
	upper_voltage: float | None = None,
) -> tuple[np.ndarray, np.ndarray, float, float, float, float]:
	mask = np.ones_like(current, dtype=bool)
	if lower_current is not None:
		mask &= current >= lower_current
	if upper_current is not None:
		mask &= current <= upper_current
	if lower_voltage is not None:
		mask &= voltage >= lower_voltage
	if upper_voltage is not None:
		mask &= voltage <= upper_voltage

	selected_current = current[mask]
	selected_voltage = voltage[mask]
	if selected_current.size < 2:
		raise ValueError("Fit region contains fewer than two points.")

	slope, intercept, slope_error, intercept_error = fit_linear(selected_current, selected_voltage)
	return selected_current, selected_voltage, slope, intercept, slope_error, intercept_error


def build_regions(current: np.ndarray, voltage: np.ndarray) -> dict[str, dict[str, object]]:
	positive = current > 0
	negative = current < 0

	peak_voltage = float(np.max(np.abs(voltage)))
	linear_lower = 0.11 * peak_voltage
	linear_upper = 0.70 * peak_voltage
	gap_threshold = 0.85 * peak_voltage

	positive_linear = fit_region(
		current[positive],
		voltage[positive],
		lower_voltage=linear_lower,
		upper_voltage=linear_upper,
	)
	positive_gap = fit_region(
		current[positive],
		voltage[positive],
		lower_voltage=gap_threshold,
	)

	negative_linear = fit_region(
		current[negative],
		voltage[negative],
		lower_voltage=-linear_upper,
		upper_voltage=-linear_lower,
	)
	negative_gap = fit_region(
		current[negative],
		voltage[negative],
		upper_voltage=-gap_threshold,
	)

	return {
		"positive_linear": {
			"current": positive_linear[0],
			"voltage": positive_linear[1],
			"slope": positive_linear[2],
			"intercept": positive_linear[3],
			"slope_error": positive_linear[4],
			"intercept_error": positive_linear[5],
		},
		"positive_gap": {
			"current": positive_gap[0],
			"voltage": positive_gap[1],
			"slope": positive_gap[2],
			"intercept": positive_gap[3],
			"slope_error": positive_gap[4],
			"intercept_error": positive_gap[5],
		},
		"negative_linear": {
			"current": negative_linear[0],
			"voltage": negative_linear[1],
			"slope": negative_linear[2],
			"intercept": negative_linear[3],
			"slope_error": negative_linear[4],
			"intercept_error": negative_linear[5],
		},
		"negative_gap": {
			"current": negative_gap[0],
			"voltage": negative_gap[1],
			"slope": negative_gap[2],
			"intercept": negative_gap[3],
			"slope_error": negative_gap[4],
			"intercept_error": negative_gap[5],
		},
		"gap_threshold": gap_threshold,
		"linear_lower": linear_lower,
		"linear_upper": linear_upper,
	}


def main() -> None:
	current, voltage = load_data(DATA_FILE)
	regions = build_regions(current, voltage)

	plt.figure(figsize=(8, 5))
	plt.plot(current * 1e3, voltage * 1e3, color="0.7", lw=1.0, label="Messdaten")

	colors = {
		"positive_linear": "tab:blue",
		"positive_gap": "tab:orange",
		"negative_linear": "tab:green",
		"negative_gap": "tab:red",
	}
	labels = {
		"positive_linear": "positiv, linear zwischen Ic und Lücke",
		"positive_gap": "positiv, oberhalb der Lücke",
		"negative_linear": "negativ, linear zwischen Ic und Lücke",
		"negative_gap": "negativ, unterhalb der Lücke",
	}

	for name in ("positive_linear", "positive_gap", "negative_linear", "negative_gap"):
		region = regions[name]
		region_current = region["current"]
		region_voltage = region["voltage"]
		slope = float(region["slope"])
		intercept = float(region["intercept"])
		sorted_indices = np.argsort(region_current)
		fitted_current = region_current[sorted_indices]
		fitted_voltage = linear_curve(fitted_current, slope, intercept)
		plt.scatter(region_current * 1e3, region_voltage * 1e3, s=8, color=colors[name], alpha=0.8)
		plt.plot(fitted_current * 1e3, fitted_voltage * 1e3, color=colors[name], lw=2.0, label=labels[name])

	plt.axhline(regions["gap_threshold"] * 1e3, color="tab:orange", ls="--", lw=0.8, alpha=0.5)
	plt.axhline(-regions["gap_threshold"] * 1e3, color="tab:red", ls="--", lw=0.8, alpha=0.5)
	plt.xlabel("Strom $I$ in mA")
	plt.ylabel("Spannung $U$ in mV")
	#plt.title("Lineare Fits für die Shunt-Auswertung bei 4.2 K")
	plt.legend(fontsize=8)
	plt.tight_layout()

	OUTPUT_PLOT.parent.mkdir(parents=True, exist_ok=True)
	plt.savefig(OUTPUT_PLOT, dpi=300)

	print(f"Daten: {DATA_FILE}")
	print(f"Plot:  {OUTPUT_PLOT}")
	print()
	for name in ("positive_linear", "positive_gap", "negative_linear", "negative_gap"):
		slope = float(regions[name]["slope"])
		intercept = float(regions[name]["intercept"])
		slope_error = float(regions[name]["slope_error"])
		intercept_error = float(regions[name]["intercept_error"])
		print(
			f"{name}: slope = {slope:.6e} ± {slope_error:.2e} V/A, "
			f"intercept = {intercept:.6e} ± {intercept_error:.2e} V"
		)


if __name__ == "__main__":
	main()
