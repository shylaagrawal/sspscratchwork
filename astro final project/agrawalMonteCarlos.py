import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import odlib

def main():
    np.random.seed(1982)

    N_samples = 1000

    obs1_coefficients = odlib.compute_coefficients("reference_stars_1.csv")
    obs2_coefficients = odlib.compute_coefficients("reference_stars_2.csv")
    obs3_coefficients = odlib.compute_coefficients("reference_stars_3.csv")

    obs1_star_positions = pd.read_csv("reference_stars_1.csv")
    obs2_star_positions = pd.read_csv("reference_stars_2.csv")
    obs3_star_positions = pd.read_csv("reference_stars_3.csv")

    sigma_ra1, sigma_dec1 = odlib.compute_sigmas(obs1_coefficients, obs1_star_positions)
    sigma_ra2, sigma_dec2 = odlib.compute_sigmas(obs2_coefficients, obs2_star_positions)
    sigma_ra3, sigma_dec3 = odlib.compute_sigmas(obs3_coefficients, obs3_star_positions)

    data = pd.read_csv("mog_test_cases.csv")

    t1 = data["t1_jd"].iloc[0]
    t2 = data["t2_jd"].iloc[0]
    t3 = data["t3_jd"].iloc[0]
    
    ra1_deg, dec1_deg = odlib.radec_to_decimal(data["ra1"].iloc[0], data["dec1"].iloc[0])
    ra2_deg, dec2_deg = odlib.radec_to_decimal(data["ra2"].iloc[0], data["dec2"].iloc[0])
    ra3_deg, dec3_deg = odlib.radec_to_decimal(data["ra3"].iloc[0], data["dec3"].iloc[0])
    
    R1 = np.array([data["sun_x_1"].iloc[0], data["sun_y_1"].iloc[0], data["sun_z_1"].iloc[0]])
    R2 = np.array([data["sun_x_2"].iloc[0], data["sun_y_2"].iloc[0], data["sun_z_2"].iloc[0]])
    R3 = np.array([data["sun_x_3"].iloc[0], data["sun_y_3"].iloc[0], data["sun_z_3"].iloc[0]])

    ra1_draws = np.random.normal(ra1_deg, sigma_ra1, N_samples)
    dec1_draws = np.random.normal(dec1_deg, sigma_dec1, N_samples)
    
    ra2_draws = np.random.normal(ra2_deg, sigma_ra2, N_samples)
    dec2_draws = np.random.normal(dec2_deg, sigma_dec2, N_samples)
    
    ra3_draws = np.random.normal(ra3_deg, sigma_ra3, N_samples)
    dec3_draws = np.random.normal(dec3_deg, sigma_dec3, N_samples)

    orbital_elements_list = []

    for i in range(N_samples):
        try:
            p_ra1, p_dec1 = odlib.decimal_to_radec(ra1_draws[i], dec1_draws[i])
            p_ra2, p_dec2 = odlib.decimal_to_radec(ra2_draws[i], dec2_draws[i])
            p_ra3, p_dec3 = odlib.decimal_to_radec(ra3_draws[i], dec3_draws[i])

            r2, v2 = odlib.gauss_method(
                t1, t2, t3,
                p_ra1, p_dec1,
                p_ra2, p_dec2,
                p_ra3, p_dec3,
                R1, R2, R3
            )
            
            orbital_elements_list.append(odlib.orbital_elements(r2, v2))

        except Exception:
            continue

    if len(orbital_elements_list) == 0:
        print("Error: All Monte Carlo trials failed.")
        return

    orbital_elements_array = np.array(orbital_elements_list)

    means = np.mean(orbital_elements_array, axis=0)
    sigmas = np.std(orbital_elements_array, axis=0, ddof=1)

    print("\nMonte Carlo Results:")
    print("----------------------------------------")
    print(f"a (AU)      : {means[0]:12.6f} +/- {sigmas[0]:10.6f}")
    print(f"e           : {means[1]:12.6f} +/- {sigmas[1]:10.6f}")
    print(f"i (deg)     : {means[2]:12.6f} +/- {sigmas[2]:10.6f}")
    print(f"Omega (deg) : {means[3]:12.6f} +/- {sigmas[3]:10.6f}")
    print(f"omega (deg) : {means[4]:12.6f} +/- {sigmas[4]:10.6f}")
    print(f"M (deg)     : {means[5]:12.6f} +/- {sigmas[5]:10.6f}")
    print("----------------------------------------")
