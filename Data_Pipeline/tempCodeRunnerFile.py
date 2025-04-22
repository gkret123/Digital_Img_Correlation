    # Plot stacked subplots for axial strain, lateral strain, and shear strain across entire beam
    fig, axs = plt.subplots(3, 1, figsize=(20, 15))
    
    cp3 = axs[0].contourf(X_entire, Y_entire, epsilon_x_entire, 500, cmap='gist_rainbow')
    fig.colorbar(cp3, ax=axs[0], label='Axial Strain, εₓ')
    axs[0].set_xlabel("Beam length, x (m)")
    axs[0].set_ylabel("Vertical position, y (m)")
    axs[0].set_title("Axial Strain Distribution Across Entire Beam")
    
    cp4 = axs[1].contourf(X_entire, Y_entire, epsilon_y_entire, 500, cmap='gist_rainbow')
    fig.colorbar(cp4, ax=axs[1], label='Lateral Strain, εᵧ')
    axs[1].set_xlabel("Beam length, x (m)")
    axs[1].set_ylabel("Vertical position, y (m)")
    axs[1].set_title("Lateral Strain Distribution Across Entire Beam")
    
    cp5 = axs[2].contourf(X_entire, Y_entire, gamma_entire, 500, cmap='gist_rainbow')
    fig.colorbar(cp5, ax=axs[2], label='Shear Strain, γ')
    axs[2].set_xlabel("Beam length, x (m)")
    axs[2].set_ylabel("Vertical position, y (m)")
    axs[2].set_title("Shear Strain Distribution Across Entire Beam")
    
    plt.tight_layout()
    plt.show()
