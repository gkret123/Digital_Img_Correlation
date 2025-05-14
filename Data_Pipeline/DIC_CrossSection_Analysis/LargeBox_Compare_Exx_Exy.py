import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import griddata

plt.rcParams.update({'font.size': 10, 'font.family': 'sans-serif'})

def load_strain_slice(file_label, strain_type, target_x):
    file_path = f"/Users/jonasmargono/Documents/DIC Data/DIC_{file_label}.xlsx"
    df = pd.read_excel(file_path)
    df.columns = df.columns.str.strip().str.replace('"', '')
    df = df[~((df['X'] == 35.4785) & (df['Y'] == 17.7537) & (df['Z'] == 483.516))]
    df_strain = df[df[strain_type] != 0].dropna(subset=['X', 'Y', strain_type])
    
    xi = np.linspace(df_strain['X'].min(), df_strain['X'].max(), 300)
    yi = np.linspace(df_strain['Y'].min(), df_strain['Y'].max(), 300)
    Xg, Yg = np.meshgrid(xi, yi)
    Z = griddata((df_strain['X'], df_strain['Y']), df_strain[strain_type], (Xg, Yg), method='linear')
    df_long = pd.DataFrame(Z, index=yi, columns=xi).stack().reset_index()
    df_long.columns = ['Y', 'X', strain_type]
    
    df_long = df_long.dropna()
    nearest_x = df_long['X'].iloc[(np.abs(df_long['X'] - target_x)).argmin()]
    return df_long[df_long['X'] == nearest_x], nearest_x

def plot_strain_comparison(label, strain_type, poly_order, target_x, ylabel, residual_color):
    slice_1000N, x_1000N = load_strain_slice(f"{label}_1000N", strain_type, target_x)
    slice_4000N, x_4000N = load_strain_slice(f"{label}_4000N", strain_type, target_x)

    coeffs_1000N = np.polyfit(slice_1000N['Y'], slice_1000N[strain_type], poly_order)
    coeffs_4000N = np.polyfit(slice_4000N['Y'], slice_4000N[strain_type], poly_order)
    fit_1000N = np.polyval(coeffs_1000N, slice_1000N['Y'])
    fit_4000N = np.polyval(coeffs_4000N, slice_4000N['Y'])

    residuals_1000N = slice_1000N[strain_type] - fit_1000N
    residuals_4000N = slice_4000N[strain_type] - fit_4000N

    def poly_to_string(coeffs):
        terms = [f"{coeff:.2e}·y^{poly_order - i}" if poly_order - i > 1
                 else f"{coeff:.2e}·y" if poly_order - i == 1
                 else f"{coeff:.2e}" for i, coeff in enumerate(coeffs)]
        return " + ".join(terms)

    eq_1000N = poly_to_string(coeffs_1000N)
    eq_4000N = poly_to_string(coeffs_4000N)

    # --- FIT PLOT ---
    plt.figure(figsize=(8, 5))
    plt.scatter(slice_1000N[strain_type], slice_1000N['Y'], s=2, alpha=0.5, color='red', label='1000N Raw')
    plt.plot(fit_1000N, slice_1000N['Y'], linestyle='--', color='maroon', linewidth=2, label='1000N Fit')
    
    plt.scatter(slice_4000N[strain_type], slice_4000N['Y'], s=2, alpha=0.5, color='blue', label='4000N Raw')
    plt.plot(fit_4000N, slice_4000N['Y'], linestyle='--', color='navy', linewidth=2, label='4000N Fit')

    plt.text(1.02, 0.95, f"1000N Fit:\n{eq_1000N}",
             transform=plt.gca().transAxes, fontsize=8, color='maroon',
             verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    plt.text(1.02, 0.65, f"4000N Fit:\n{eq_4000N}",
             transform=plt.gca().transAxes, fontsize=8, color='navy',
             verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    plt.xlabel(f"{ylabel} (×10⁻⁴)")
    plt.ylabel("y (mm)")
    plt.title(f"{ylabel} vs y at x ≈ {x_1000N:.2f} mm (Combined {title} Tests)")
    plt.grid(True)
    plt.legend(markerscale=2, loc='lower right', fontsize=9)
    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x*1e4:.2f}"))
    plt.tight_layout()
    plt.xlim(left=min(slice_1000N[strain_type].min(), slice_4000N[strain_type].min()) * 0.95,
             right=max(slice_1000N[strain_type].max(), slice_4000N[strain_type].max()) * 1.05)
    plt.ylim(slice_1000N['Y'].min(), slice_1000N['Y'].max())
    plt.show()

    # --- RESIDUAL PLOT ---
    plt.figure(figsize=(6, 3))
    plt.scatter(residuals_1000N, slice_1000N['Y'], s=5, alpha=0.6, color=residual_color, label='1000N Residual')
    plt.scatter(residuals_4000N, slice_4000N['Y'], s=5, alpha=0.6, color='blue', label='4000N Residual')
    plt.axvline(0, color='black', linestyle='--')
    plt.xlabel(f"Residual {ylabel} (×10⁻⁴)")
    plt.ylabel("y (mm)")
    plt.title(f"Residuals of {ylabel} Fit at x ≈ {x_1000N:.2f} mm ({title})")
    plt.grid(True)
    plt.legend(fontsize=9)
    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x*1e4:.2f}"))
    plt.tight_layout()
    plt.show()

# Example usage
title = "Large Beam"
label = "LargeBox"
plot_strain_comparison(label, strain_type="exx", poly_order=1, target_x=30, ylabel="εₓₓ", residual_color='red')
#plot_strain_comparison(label, strain_type="exy", poly_order=2, target_x=30, ylabel="εₓᵧ", residual_color='red')
