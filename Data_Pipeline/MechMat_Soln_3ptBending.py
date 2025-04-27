"""
3-Point Bending Beam Analysis

Author: Gabriel Kret 
Instructors: Dr. David Wootton and Dr. Kamau Wright
Course: ME-360-A: Experimentation
Date: Spring 2025


This script analyzes a simply supported beam subjected to a central load using Mechanics of Materials beam theory.
It computes the bending moment, shear force, deflection, and strain/stress distributions across the beam's cross-section.
It also visualizes the results using plots and heatmaps.
This script is intended to be used alongside the Digital Image Correlation (DIC) analysis of a 3-point bending test.
The results from this script can be compared with the DIC results to validate the mechanical behavior of the beam.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

"""
F_applied and X_coord are the user inputs for the loading condition at a specific x coordinate along the beam.
enter a value for F_applied to specify the loading condition you wish to analyze.
The value of F_applied is the force applied at the center of the beam.
The script will use this value to calculate the bending moment, shear forces, deflections, etc.
NOTE: It will also be used to approximate whether the case puts the beam in the elastic or plastic region.

enter a value for X_coord as the position along the beam where you want to analyze the bending moment, shear force, and deflection.
The value of X_coord is the distance from the left end of the beam to the point where you want to analyze the bending moment, shear force, and deflection.
"""
F_applied = 4000  # Applied force in Newtons (example value)
x_coord = 0.03 # Position along the beam in meters (example value, 6.25 cm)

# Beam geometry and material properties

L = 0.12  # Beam length in meters (12 cm) (measured)
width = 0.019  # meters (measured)
height = 0.019 # meters (measured)
        
# Note: The inner dimensions are only used if the beam is hollow. write 'None' if the beam is solid.
inner_width = 0.013 #measured
inner_height = inner_width  # meters 

#known/estimated material properties
poisson = 0.33  # Poisson's ratio (material property)
E = 69e9      # Young's modulus in Pascals (example value for Aluminum 6061) (Material property)
yield_strength = 276e6  # Yield strength in Pascals (example value for Aluminum 6061) (Material property)
    
G = E / (2 * (1 + poisson)) # Compute the shear modulus from Young's modulus and Poisson's ratio

def compute_inertia(width, height, inner_width = None, inner_height = None):
    """
    Compute the second moment of area (I) for a rectangular or hollow rectangular cross-section.
    """
    if inner_width is not None and inner_height is not None:
        # Hollow section
        I_outer = (width * height**3) / 12.0
        I_inner = (inner_width * inner_height**3) / 12.0
        I = I_outer - I_inner
    else:
        # Solid section
        I = (width * height**3) / 12.0
    return I
   

def bending_moment(x, L, F_applied):
    """
    Compute the bending moment at a distance x along a simply supported beam 
    with a central load F_applied.
    For x <= L/2: M = F_applied*x/2; for x > L/2: M = F_applied*(L-x)/2.
    """
    if x <= L / 2:
        return F_applied * x / 2.0
    else:
        return F_applied * (L - x) / 2.0

def bending_moment_array(x_array, L, F_applied):
    """ Vectorized bending moment over an array of x coordinates """
    return np.array([bending_moment(x, L, F_applied) for x in x_array])

def shear_force(x, L, F_applied):
    """
    Compute the shear force at a distance x along the beam.
    For a simply supported beam with a central load:
      if x < L/2, shear V = F_applied/2; if x > L/2, V = -F_applied/2.
    (A jump exists at the midspan.)
    """
    if x < L / 2:
        return F_applied / 2.0
    else:
        return -F_applied / 2.0

def shear_force_array(x_array, L, F_applied):
    """  shear force over an array of x coordinates """
    return np.array([shear_force(x, L, F_applied) for x in x_array])

def deflection(x, L, F_applied, E, I):
    """
    Compute the vertical deflection at position x along the beam.
    For a simply supported beam with a central load, the deflection is given by:
      For 0 <= x <= L/2:
         v(x) = (F_applied*x*(3*L**2 - 4*x**2))/(48*E*I)
      For L/2 <= x <= L:
         v(x) = (F_applied*(L-x)*(3*L**2 - 4*(L-x)**2))/(48*E*I)
    """
    if x <= L / 2:
        return (F_applied * x * (3 * L**2 - 4 * x**2)) / (48 * E * I)
    else:
        return (F_applied * (L - x) * (3 * L**2 - 4 * (L - x)**2)) / (48 * E * I)

def deflection_array(x_array, L, F_applied, E, I):
    """ Vectorized deflection calculation along the beam """
    return np.array([deflection(x, L, F_applied, E, I) for x in x_array])

def strain_distribution(y, M, E, I):
    """
    Compute the axial strain at a point y (distance from the neutral axis) due to bending.
    εₓ(y) = - y * M / (E * I)
    """
    return - y * M / (E * I)

def shear_stress_distribution(y, V, width, height):
    """
    Compute the shear stress distribution in a rectangular cross-section at the neutral axis.
    For a rectangular beam:
      τ(y) = (3/2)*(V/(width*height))*(1 - (2*y/height)**2)
    """
    return (3/2) * (V / (width * height)) * (1 - (2 * y / height)**2)
    
def elastic_or_plastic(F_applied, I):
    
    #Determine if the beam is in the elastic or plastic region for the force case applied

    c = height / 2  # Distance from the neutral axis to the outer fiber
    # Calculate the maximum bending stress allowed
    M_yield = yield_strength * I / c  # Yield moment
    # Calculate the maximum load that can be applied without yielding
    P_yield = 4 * M_yield / L  # Yield load, for a simply supported beam with a central load
    # Calculate the maximum bending stress induced by the applied load
    M = bending_moment(L / 2, L, F_applied)  # Maximum moment at the center of the beam
    sigma_max = M * c / I  # Maximum bending stress induced by the applied load

    print("The critical load is: ", P_yield, "N")
    print("The applied load is: ", F_applied, "N")
    print("The maximum bending stress from the applied load is: ", sigma_max/10**6, "MPa")
    print("The yield strength is: ", yield_strength/10**6, "MPa")

    #compare with yield strength
    if sigma_max > yield_strength:
        print("\nThe beam is in the plastic region, the force applied induced a bending stress higher than the yield strength. \n")

        print(f"The maximum bending stress is: {sigma_max/10**6} MPa, which is greater than the yield strength of {yield_strength/10**6} MPa.")
        print(f"The applied load of {F_applied} N is greater than the critical load of {P_yield} N, The beam is expected to yield at an applied load of {P_yield} N.")
        return "Plastic"
    else:
        print("\nThe beam is in the elastic region, the force applied does not induce a bending stress higher than the yield strength. \n")

        print(f"The maximum bending stress is: {sigma_max/10**6} MPa, which is less than the yield strength of {yield_strength/10**6} MPa.")
        print(f"An applied load of {P_yield} N would cause the beam to cross into the plastic region.")
        return "Elastic"
    


def main():
    # the user inputs for the loading condition at a specific x coordinate along the beam (Force and x position) are defined above
    
    
    # Calculate the moment of inertia
    I = compute_inertia(width, height, inner_width, inner_height)
    # Bending moment at the chosen x coordinate
    M = bending_moment(x_coord, L, F_applied)
    # Shear force at the chosen x coordinate
    V = shear_force(x_coord, L, F_applied)
    
    
    
    # Cross-Section Analysis
    
    # Define vertical positions (y) across the cross-section (neutral axis at y = 0)
    y = np.linspace(-height / 2, height / 2, 100)
    # Compute axial strain distribution at the given x coordinate
    epsilon_x = strain_distribution(y, M, E, I)
    # Lateral strain (due to Poisson effect)
    epsilon_y = - poisson * epsilon_x
    # Bending stress distribution (stress = E * strain)
    sigma = E * epsilon_x
    # Shear stress distribution across the height
    tau = shear_stress_distribution(y, V, width, height)
    
    print("\nSample strain/stress values (every 10th point) at x = {:.4f} m:".format(x_coord))
    for i in range(0, len(y), 10):
        print(f"y = {y[i]:.4f} m, εₓ = {epsilon_x[i]:.6e}, εᵧ = {epsilon_y[i]:.6e}, σ = {sigma[i]:.6e} Pa, τ = {tau[i]:.6e} Pa")
    
    print("\n--- Beam Analysis Results ---\n")
    #bending region
    print("Understanding the beams elastic and plastic regions:\n")
    elastic_or_plastic(F_applied, I)

    print("\nBeam Geometry and Material Properties:\n")

    print(f"Moment of Inertia (I): {I:.6e} m^4")
    print(f"Bending Moment (M) at x = {x_coord:.4f} m: {M:.6e} Nm")
    print(f"Shear Force (V) at x = {x_coord:.4f} m: {V:.6e} N")
    
    # Maximum deflection for a simply supported beam with a central load
    delta_max = F_applied * L**3 / (48 * E * I)
    print(f"Maximum Deflection: {delta_max:.6e} m")
    


    # Plot 1–3: Strains vs. y in cross‐section (axial, lateral, shear)
    fig, axs = plt.subplots(1, 3, figsize=(18, 4))
    # Axial Strain
    axs[0].plot(epsilon_x, y, 'b-', label='Axial Strain, εₓ')
    axs[0].set_xlabel("Axial Strain, εₓ")
    axs[0].set_ylabel("Vertical position, y (m)")
    axs[0].set_title(f"Axial Strain Distribution\n(at x = {x_coord:.4f} m)")
    axs[0].grid(True)
    axs[0].legend()

    # Lateral (bending) Strain
    axs[1].plot(epsilon_y, y, 'g-', label='Lateral Strain, εᵧ')
    axs[1].set_xlabel("Lateral Strain, εᵧ")
    axs[1].set_ylabel("Vertical position, y (m)")
    axs[1].set_title(f"Lateral Strain Distribution\n(at x = {x_coord:.4f} m)")
    axs[1].grid(True)
    axs[1].legend()

    # Shear Strain
    gamma = tau / G
    axs[2].plot(gamma, y, 'm-', label='Shear Strain, εₓᵧ')
    axs[2].set_xlabel("Shear Strain, εₓᵧ")
    axs[2].set_ylabel("Vertical position, y (m)")
    axs[2].set_title(f"Shear Strain Distribution\n(at x = {x_coord:.4f} m)")
    axs[2].grid(True)
    axs[2].legend()

    plt.tight_layout()
    plt.suptitle(f"Strain Distributions in Cross-Section\nat x = {x_coord:.4f} m", fontsize=16)
    
    # Plot 4-6: Stress distributions vs. y in cross‐section
    fig, axs = plt.subplots(1, 3, figsize=(18, 4))

    # Normal bending stress 
    axs[0].plot(sigma, y, 'r-', label='σₓ (Pa)')
    axs[0].set_xlabel("Normal Stress, σₓ (Pa)")
    axs[0].set_ylabel("Vertical position, y (m)")
    axs[0].set_title(f"Normal Stress Distribution\n(at x = {x_coord:.4f} m)")
    axs[0].grid(True)
    axs[0].legend()

    # Lateral stress (Poisson effect)
    sigma_y = E * epsilon_y
    axs[1].plot(sigma_y, y, 'g-', label='σᵧ (Pa)')
    axs[1].set_xlabel("Lateral Stress, σᵧ (Pa)")
    axs[1].set_ylabel("Vertical position, y (m)")
    axs[1].set_title(f"Lateral Stress Distribution\n(at x = {x_coord:.4f} m)")
    axs[1].grid(True)
    axs[1].legend()

    # Shear stress 
    axs[2].plot(tau, y, 'm-', label='τ (Pa)')
    axs[2].set_xlabel("Shear Stress, τ (Pa)")
    axs[2].set_ylabel("Vertical position, y (m)")
    axs[2].set_title(f"Shear Stress Distribution\n(at x = {x_coord:.4f} m)")
    axs[2].grid(True)
    axs[2].legend()

    plt.tight_layout()
    plt.suptitle(f"Stress Distributions in Cross-Section\nat x = {x_coord:.4f} m", fontsize=16)
    
    # Plots 6-9: Heatmaps at the cross‐section specified
    ny_cs, nz = 50, 50
    y_cs = np.linspace(-height/2, height/2, ny_cs)
    z_cs = np.linspace(-width/2, width/2, nz)
    Z_cs, Y_cs = np.meshgrid(z_cs, y_cs)

    # Strain fields
    exx_cs = strain_distribution(Y_cs, M, E, I)
    eyy_cs = - poisson * exx_cs

    # Shear strain 
    V_cs = shear_force(x_coord, L, F_applied)
    tau_cs = shear_stress_distribution(Y_cs, V_cs, width, height)
    gamma_cs = tau_cs / G
    epsilon_x_y = gamma_cs / 2
    fig, axs = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)

    # Axial Strain 
    cp0 = axs[0].contourf(Z_cs, Y_cs, exx_cs, 200, cmap='gist_rainbow')
    fig.colorbar(cp0, ax=axs[0], label='Axial Strain, εₓ')
    axs[0].set_xlabel("Horizontal coordinate, z (m)")
    axs[0].set_ylabel("Vertical coordinate, y (m)")
    axs[0].set_title(f"Axial Strain (εₓ)\nCross-Section at x = {x_coord:.4f} m")

    # Lateral Strain
    cp1 = axs[1].contourf(Z_cs, Y_cs, eyy_cs, 200, cmap='gist_rainbow')
    fig.colorbar(cp1, ax=axs[1], label='Lateral Strain, εᵧ')
    axs[1].set_xlabel("Horizontal coordinate, z (m)")
    axs[1].set_ylabel("Vertical coordinate, y (m)")
    axs[1].set_title(f"Lateral Strain (εᵧ)\nCross-Section at x = {x_coord:.4f} m")

    # Shear Strain 
    #cp2 = axs[2].contourf(Z_cs, Y_cs, gamma_cs, 200, cmap='gist_rainbow')
    cp2 = axs[2].contourf(Z_cs, Y_cs, epsilon_x_y, 200, cmap='gist_rainbow') #changed to include exy to match VIC data
    fig.colorbar(cp2, ax=axs[2], label='Shear Strain, εₓᵧ')
    axs[2].set_xlabel("Horizontal coordinate, z (m)")
    axs[2].set_ylabel("Vertical coordinate, y (m)")
    axs[2].set_title(f"Shear Strain (εₓᵧ)\nCross-Section at x = {x_coord:.4f} m")
    
    plt.suptitle(f"Heatmap Display of Strain Distributions in Cross-Section\nat x = {x_coord:.4f} m", fontsize=16)


    # Beam Diagram Plots (along x)
    x_beam = np.linspace(0, L, 200)
    M_beam = bending_moment_array(x_beam, L, F_applied)
    V_beam = shear_force_array(x_beam, L, F_applied)
    v_beam = deflection_array(x_beam, L, F_applied, E, I)
    
    
    # Plot 10-11: Shear and bending moment diagrams
    fig, ax = plt.subplots(2, figsize=(6,4))
    ax[0].plot(x_beam, V_beam, 'c-', label='Shear Force, V(x)')
    ax[0].set_xlabel("Beam length, x (m)")
    ax[0].set_ylabel("Shear Force, V (N)")
    ax[0].set_title("Shear Force Diagram")
    ax[0].grid(True)
    ax[0].legend()
    ax[1].plot(x_beam, M_beam, 'k-', label='Bending Moment, M(x)')
    ax[1].set_xlabel("Beam length, x (m)")
    ax[1].set_ylabel("Bending Moment, M (Nm)")
    ax[1].set_title("Bending Moment Diagram")
    ax[1].grid(True)
    ax[1].legend()
    plt.tight_layout()
    
    
    # Plot 12: Deflection Curve
    plt.figure(figsize=(12,10))
    plt.plot(x_beam, v_beam, 'g-', label='Deflection, v(x)')
    plt.xlabel("Beam length, x (m)")
    plt.ylabel("Vertical Deflection, v (m)")
    plt.title("Deflection Curve of the Beam")
    plt.grid(True)
    plt.legend()
    
    
    # Plots 13-15: stacked heatmap subplots for axial strain, lateral strain, and shear strain across entire beam.
    # Create a grid covering the entire beam: x from 0 to L and y from -height/2 to height/2.
    nx_entire, ny_entire = 2000, 500
    x_entire = np.linspace(0, L, nx_entire)
    y_entire = np.linspace(-height/2, height/2, ny_entire)
    X_entire, Y_entire = np.meshgrid(x_entire, y_entire)
    # Compute bending moment at each x location (note: M is independent of y)
    M_entire = np.array([bending_moment(x, L, F_applied) for x in x_entire])
    # Broadcast M_entire along y-direction
    M_grid = np.tile(M_entire, (ny_entire, 1))
    # Compute axial strain distribution over the entire beam:
    epsilon_x_entire = - Y_entire * M_grid / (E * I)
    # Lateral strain distribution using Poisson's ratio:
    epsilon_y_entire = - poisson * epsilon_x_entire
    # shear strain distribution using shear force:

    V_entire = np.array([shear_force(x, L, F_applied) for x in x_entire])
    V_grid = np.tile(V_entire, (ny_entire, 1))
    # First compute shear stress distribution
    tau_entire = shear_stress_distribution(Y_entire, V_grid, width, height)
    # Convert shear stress to shear strain
    gamma_entire = tau_entire / G
    epsilon_x_y_entire = gamma_entire / 2  # Shear strain (engineering shear strain)

    # create subplots with automatic spacing for titles/colorbars
    fig, axs = plt.subplots(3, 1, figsize=(16, 10), constrained_layout=True)

    # Axial strain
    cp3 = axs[0].contourf(X_entire, Y_entire, epsilon_x_entire, 500, cmap='gist_rainbow')
    cb0 = fig.colorbar(cp3, ax=axs[0], label='Axial Strain, εₓ')
    vmin0, vmax0 = epsilon_x_entire.min(), epsilon_x_entire.max()
    ticks0 = np.linspace(vmin0, vmax0, 8)
    cb0.set_ticks([vmin0, *ticks0, vmax0])
    cb0.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
    axs[0].set_xlabel("Beam length, x (m)")
    axs[0].set_ylabel("Vertical position, y (m)")
    axs[0].set_title("Axial Strain Distribution Across Entire Beam", pad=12)

    # Lateral strain
    cp4 = axs[1].contourf(X_entire, Y_entire, epsilon_y_entire, 500, cmap='gist_rainbow')
    cb1 = fig.colorbar(cp4, ax=axs[1], label='Lateral Strain, εᵧ')
    vmin1, vmax1 = epsilon_y_entire.min(), epsilon_y_entire.max()
    ticks1 = np.linspace(vmin1, vmax1, 8)
    cb1.set_ticks([vmin1, *ticks1, vmax1])
    cb1.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
    axs[1].set_xlabel("Beam length, x (m)")
    axs[1].set_ylabel("Vertical position, y (m)")
    axs[1].set_title("Lateral Strain Distribution Across Entire Beam", pad=12)

    # Shear strain
    #cp5 = axs[2].contourf(X_entire, Y_entire, gamma_entire, 500, cmap='gist_rainbow')
    cp5 = axs[2].contourf(X_entire, Y_entire, epsilon_x_y_entire, 500, cmap='gist_rainbow') #changed to include exy to match VIC data
    cb2 = fig.colorbar(cp5, ax=axs[2], label='Shear Strain, εₓᵧ')
    #vmin2, vmax2 = gamma_entire.min(), gamma_entire.max()
    vmin2, vmax2 = epsilon_x_y_entire.min(), epsilon_x_y_entire.max() #changed to include exy to match VIC data
    ticks2 = np.linspace(vmin2, vmax2, 8)
    cb2.set_ticks([vmin2, *ticks2, vmax2])
    cb2.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
    axs[2].set_xlabel("Beam length, x (m)")
    axs[2].set_ylabel("Vertical position, y (m)")
    axs[2].set_title("Shear Strain Distribution Across Entire Beam", pad=12)

    # 5) show
    plt.show()

if __name__ == "__main__":
    main()
