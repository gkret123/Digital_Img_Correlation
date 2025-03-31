import numpy as np
import matplotlib.pyplot as plt

# Beam geometry and material properties
L = 0.12  # Beam length in meters (12 cm)
width = 1.5 / 39.3700787  # Convert 1.5 inches to meters
height = 1.5 / 39.3700787  # Convert 1.5 inches to meters
poisson = 0.33  # Poisson's ratio 
E = 69e9      # Young's modulus in Pascals (example value for Aluminum 6061)

def compute_inertia(width, height):
    """
    Compute the second moment of area (I) for a rectangular cross-section.
    I = (width * height^3) / 12.
    """
    return (width * height**3) / 12.0

def bending_moment(x, L, F):
    """
    Compute the bending moment at a distance x along a simply supported beam 
    with a central load F.
    For x <= L/2: M = F*x/2; for x > L/2: M = F*(L-x)/2.
    """
    if x <= L / 2:
        return F * x / 2.0
    else:
        return F * (L - x) / 2.0

def bending_moment_array(x_array, L, F):
    """ Vectorized bending moment over an array of x coordinates """
    return np.array([bending_moment(x, L, F) for x in x_array])

def shear_force(x, L, F):
    """
    Compute the shear force at a distance x along the beam.
    For a simply supported beam with a central load:
      if x < L/2, shear V = F/2; if x > L/2, V = -F/2.
    (A jump exists at the midspan.)
    """
    if x < L / 2:
        return F / 2.0
    else:
        return -F / 2.0

def shear_force_array(x_array, L, F):
    """ Vectorized shear force over an array of x coordinates """
    return np.array([shear_force(x, L, F) for x in x_array])

def deflection(x, L, F, E, I):
    """
    Compute the vertical deflection at position x along the beam.
    For a simply supported beam with a central load, the deflection is given by:
      For 0 <= x <= L/2:
         v(x) = (F*x*(3*L**2 - 4*x**2))/(48*E*I)
      For L/2 <= x <= L:
         v(x) = (F*(L-x)*(3*L**2 - 4*(L-x)**2))/(48*E*I)
    """
    if x <= L / 2:
        return (F * x * (3 * L**2 - 4 * x**2)) / (48 * E * I)
    else:
        return (F * (L - x) * (3 * L**2 - 4 * (L - x)**2)) / (48 * E * I)

def deflection_array(x_array, L, F, E, I):
    """ Vectorized deflection calculation along the beam """
    return np.array([deflection(x, L, F, E, I) for x in x_array])

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

def main():
    # User inputs for the loading condition at a specific x coordinate along the beam
    F = float(input("Enter the applied force (N): "))
    x_coord = float(input("Enter the x coordinate along the beam (m): "))
    
    # Calculate the moment of inertia
    I = compute_inertia(width, height)
    # Bending moment at the chosen x coordinate
    M = bending_moment(x_coord, L, F)
    # Shear force at the chosen x coordinate
    V = shear_force(x_coord, L, F)
    
    print("\n--- Beam Analysis Results ---")
    print(f"Moment of Inertia (I): {I:.6e} m^4")
    print(f"Bending Moment (M) at x = {x_coord:.4f} m: {M:.6e} N·m")
    print(f"Shear Force (V) at x = {x_coord:.4f} m: {V:.6e} N")
    
    # Maximum deflection for a simply supported beam with a central load
    delta_max = F * L**3 / (48 * E * I)
    print(f"Maximum Deflection: {delta_max:.6e} m")
    
    # =======================
    # Cross-Section Analysis (at x = x_coord)
    # =======================
    # Define vertical positions (y) across the cross-section (neutral axis at y = 0)
    y = np.linspace(-height / 2, height / 2, 100)
    # Compute axial strain distribution at the given x coordinate
    epsilon_x = strain_distribution(y, M, E, I)
    # Lateral strain (due to Poisson effect)
    epsilon_y = - poisson * epsilon_x
    # Bending stress distribution (σ = E * εₓ)
    sigma = E * epsilon_x
    # Shear stress distribution across the height
    tau = shear_stress_distribution(y, V, width, height)
    
    # Display sample values every 5 points
    print("\nSample strain/stress values (every 5th point) at x = {:.4f} m:".format(x_coord))
    for i in range(0, len(y), 5):
        print(f"y = {y[i]:.4f} m, εₓ = {epsilon_x[i]:.6e}, εᵧ = {epsilon_y[i]:.6e}, σ = {sigma[i]:.6e} Pa, τ = {tau[i]:.6e} Pa")
    
    # ------------------------------
    # Plot 1: Axial Strain vs. y (at x = x_coord)
    plt.figure(figsize=(6,4))
    plt.plot(epsilon_x, y, 'b-', label='Axial Strain, εₓ')
    plt.xlabel("Axial Strain, εₓ")
    plt.ylabel("Vertical position, y (m)")
    plt.title("Axial Strain Distribution Through Cross-Section\n(at x = {:.4f} m)".format(x_coord))
    plt.grid(True)
    plt.legend()
    
    # ------------------------------
    # Plot 2: Bending Stress vs. y (at x = x_coord)
    plt.figure(figsize=(6,4))
    plt.plot(sigma, y, 'r-', label='Bending Stress, σ')
    plt.xlabel("Bending Stress, σ (Pa)")
    plt.ylabel("Vertical position, y (m)")
    plt.title("Bending Stress Distribution Through Cross-Section\n(at x = {:.4f} m)".format(x_coord))
    plt.grid(True)
    plt.legend()
    
    # ------------------------------
    # Plot 3: Shear Stress vs. y (at x = x_coord)
    plt.figure(figsize=(6,4))
    plt.plot(tau, y, 'm-', label='Shear Stress, τ')
    plt.xlabel("Shear Stress, τ (Pa)")
    plt.ylabel("Vertical position, y (m)")
    plt.title("Shear Stress Distribution Through Cross-Section\n(at x = {:.4f} m)".format(x_coord))
    plt.grid(True)
    plt.legend()
    
    # ------------------------------
    # Plot 4 & 5: Heatmaps of εₓ and εᵧ in the cross-section (at x = x_coord)
    # Create a grid for the cross-section (y: vertical, z: horizontal)
    ny_cs, nz = 50, 50
    y_cs = np.linspace(-height/2, height/2, ny_cs)
    z_cs = np.linspace(-width/2, width/2, nz)
    Z_cs, Y_cs = np.meshgrid(z_cs, y_cs)  # Y_cs varies vertically
    # For pure bending, εₓ is only a function of y (constant along z)
    exx_cs = strain_distribution(Y_cs, M, E, I)
    eyy_cs = - poisson * exx_cs
    
    plt.figure(figsize=(6,5))
    cp = plt.contourf(Z_cs, Y_cs, exx_cs, 200, cmap='viridis')
    plt.colorbar(cp, label='Axial Strain, εₓ')
    plt.xlabel("Horizontal coordinate, z (m)")
    plt.ylabel("Vertical coordinate, y (m)")
    plt.title("Heatmap of Axial Strain (εₓ)\nin Cross-Section (at x = {:.4f} m)".format(x_coord))
    
    plt.figure(figsize=(6,5))
    cp2 = plt.contourf(Z_cs, Y_cs, eyy_cs, 200, cmap='plasma')
    plt.colorbar(cp2, label='Lateral Strain, εᵧ')
    plt.xlabel("Horizontal coordinate, z (m)")
    plt.ylabel("Vertical coordinate, y (m)")
    plt.title("Heatmap of Lateral Strain (εᵧ)\nin Cross-Section (at x = {:.4f} m)".format(x_coord))
    
    # ===============================
    # Beam Diagram Plots (along x)
    # ===============================
    x_beam = np.linspace(0, L, 200)
    M_beam = bending_moment_array(x_beam, L, F)
    V_beam = shear_force_array(x_beam, L, F)
    v_beam = deflection_array(x_beam, L, F, E, I)
    
    # ------------------------------
    # Plot 6: Bending Moment Diagram
    plt.figure(figsize=(6,4))
    plt.plot(x_beam, M_beam, 'k-', label='Bending Moment, M(x)')
    plt.xlabel("Beam length, x (m)")
    plt.ylabel("Bending Moment, M (N·m)")
    plt.title("Bending Moment Diagram")
    plt.grid(True)
    plt.legend()
    
    # ------------------------------
    # Plot 7: Shear Force Diagram
    plt.figure(figsize=(6,4))
    plt.plot(x_beam, V_beam, 'c-', label='Shear Force, V(x)')
    plt.xlabel("Beam length, x (m)")
    plt.ylabel("Shear Force, V (N)")
    plt.title("Shear Force Diagram")
    plt.grid(True)
    plt.legend()
    
    # ------------------------------
    # Plot 8: Deflection Curve
    plt.figure(figsize=(6,4))
    plt.plot(x_beam, v_beam, 'g-', label='Deflection, v(x)')
    plt.xlabel("Beam length, x (m)")
    plt.ylabel("Vertical Deflection, v (m)")
    plt.title("Deflection Curve of the Beam")
    plt.grid(True)
    plt.legend()
    
    # ------------------------------
    # NEW: Heatmaps of Strains Across the Entire Beam
    # Create a grid covering the entire beam: x from 0 to L and y from -height/2 to height/2.
    nx_entire, ny_entire = 2000, 500
    x_entire = np.linspace(0, L, nx_entire)
    y_entire = np.linspace(-height/2, height/2, ny_entire)
    X_entire, Y_entire = np.meshgrid(x_entire, y_entire)
    # Compute bending moment at each x location (note: M is independent of y)
    M_entire = np.array([bending_moment(x, L, F) for x in x_entire])
    # Broadcast M_entire along y-direction
    M_grid = np.tile(M_entire, (ny_entire, 1))
    # Compute axial strain distribution over the entire beam:
    epsilon_x_entire = - Y_entire * M_grid / (E * I)
    # Lateral strain distribution using Poisson's ratio:
    epsilon_y_entire = - poisson * epsilon_x_entire
    
    # Plot 9: Heatmap of Axial Strain (εₓ) across the entire beam
    plt.figure(figsize=(20,5))
    cp3 = plt.contourf(X_entire, Y_entire, epsilon_x_entire, 500, cmap='viridis')
    plt.colorbar(cp3, label='Axial Strain, εₓ')
    plt.xlabel("Beam length, x (m)")
    plt.ylabel("Vertical position, y (m)")
    plt.title("Axial Strain Distribution Across Entire Beam")
    
    # Plot 10: Heatmap of Lateral Strain (εᵧ) across the entire beam
    plt.figure(figsize=(20,5))
    cp4 = plt.contourf(X_entire, Y_entire, epsilon_y_entire, 500, cmap='plasma')
    plt.colorbar(cp4, label='Lateral Strain, εᵧ')
    plt.xlabel("Beam length, x (m)")
    plt.ylabel("Vertical position, y (m)")
    plt.title("Lateral Strain Distribution Across Entire Beam")
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
