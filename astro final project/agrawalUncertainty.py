import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def main():

    # Set random seed for reproducibility
    np.random.seed(42)

    # Sample input parameters and observations
    N_samples = 1000

    # Base observations (t, alpha, delta) and their 1-sigma uncertainties
    # alpha and delta uncertainties derived from LSPR residuals
    obs_alpha = np.array([10.123, 11.456, 12.789])
    obs_delta = np.array([45.123, 46.456, 47.789])
    sigma_alpha = np.array([0.0012, 0.0015, 0.0011])
    sigma_delta = np.array([0.0014, 0.0013, 0.0016])

    # Placeholder OD function: replace with your actual OD solver call
    def run_od_solver(ra, dec):
        # Simulating non-linear responses and correlations between parameters
        a = 1.25 + 0.02 * (ra[0] - 10.123) + np.random.normal(0, 0.005)
        e = 0.15 + 0.01 * (dec[0] - 45.123) + np.random.normal(0, 0.002)
        i = 12.4 + 0.5 * (ra[1] - 11.456) + np.random.normal(0, 0.05)
        omega = 145.2 + 1.2 * (dec[1] - 46.456) + np.random.normal(0, 0.2)
        Omega = 88.3 + 0.8 * (ra[2] - 12.789) + np.random.normal(0, 0.1)
        Tp = 2460000.5 + 0.1 * (dec[2] - 47.789) + np.random.normal(0, 0.01)
        return [a, e, i, omega, Omega, Tp]

    # Step 1: Perform Monte Carlo sampling
    orbital_element_samples = []

    for _ in range(N_samples):
        # Perturb inputs using normal distribution
        perturbed_alpha = np.random.normal(obs_alpha, sigma_alpha)
        perturbed_delta = np.random.normal(obs_delta, sigma_delta)
        
        # Run OD code with perturbed observations
        elements = run_od_solver(perturbed_alpha, perturbed_delta)
        orbital_element_samples.append(elements)

    # Convert results into a numpy array of shape (N_samples, 6)
    data = np.array(orbital_element_samples)

    labels = [r'$a$ (AU)', r'$e$', r'$i$ (deg)', r'$\omega$ (deg)', r'$\Omega$ (deg)', r'$T_p$ (JD)']
    num_params = len(labels)

    # Step 2: Construct the Corner Plot grid layout
    fig = plt.figure(figsize=(12, 12))
    gs = gridspec.GridSpec(num_params, num_params, hspace=0.1, wspace=0.1)

    for row in range(num_params):
        for col in range(num_params):
            # Upper triangle cells are turned off
            if col > row:
                ax = fig.add_subplot(gs[row, col])
                ax.axis('off')
                continue
                
            ax = fig.add_subplot(gs[row, col])
            
            # Main Diagonal: 1D Histograms
            if row == col:
                ax.hist(data[:, row], bins=30, color='skyblue', edgecolor='black', histtype='stepfilled')
                ax.set_yticks([])
                if row < num_params - 1:
                    ax.set_xticklabels([])
                else:
                    ax.set_xlabel(labels[row])
                    
            # Lower Triangle: 2D Scatter plots
            else:
                ax.scatter(data[:, col], data[:, row], s=2, alpha=0.4, color='navy')
                
                # Manage label visibility to clean up shared axes
                if row == num_params - 1:
                    ax.set_xlabel(labels[col])
                else:
                    ax.set_xticklabels([])
                    
                if col == 0:
                    ax.set_ylabel(labels[row])
                else:
                    ax.set_yticklabels([])

    plt.tight_layout()
    plt.show()