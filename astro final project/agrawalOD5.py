# Shyla Agrawal
# 7/21/2026
# SSP Orbit Determination Code 5

import numpy as np
import odlib
import pandas as pd
from importlib import reload

reload(odlib)

def main():

    data = pd.read_csv("mog_test_cases.csv")

    for index, row in data.iterrows():

        print(f"Asteroid: {row['asteroid_name']}")

        t1 = row["t1_jd"]
        t2 = row["t2_jd"]
        t3 = row["t3_jd"]
        ra1 = row["ra1"]
        dec1 = row["dec1"]
        ra2 = row["ra2"]
        dec2 = row["dec2"]
        ra3 = row["ra3"]
        dec3 = row["dec3"]

        R1 = np.array([row["sun_x_1"], row["sun_y_1"], row["sun_z_1"]])
        R2 = np.array([row["sun_x_2"], row["sun_y_2"], row["sun_z_2"]])
        R3 = np.array([row["sun_x_3"], row["sun_y_3"], row["sun_z_3"]])
        R1 = odlib.ecliptic_to_equatorial(R1)
        R2 = odlib.ecliptic_to_equatorial(R2)
        R3 = odlib.ecliptic_to_equatorial(R3)

        r2, v2 = odlib.gauss_method(
            t1, t2, t3,
            ra1, dec1,
            ra2, dec2,
            ra3, dec3,
            R1, R2, R3,
            order_four=True,
            tolerance=1e-10
        )
        
        # Orbital elements
        r2 = odlib.equatorial_to_ecliptic(r2)
        v2 = odlib.equatorial_to_ecliptic(v2)

        a, e, i, Omega, omega, M = odlib.orbital_elements(r2, v2)


        jpl_r2 = np.array([float(row["jpl_r2_x"]), float(row["jpl_r2_y"]), float(row["jpl_r2_z"])])
        jpl_v2 = np.array([ row["jpl_v2_x"], row["jpl_v2_y"],row["jpl_v2_z"]])

        print("\nOrbital Element Error")
        jpl_a = row["jpl_a"]
        jpl_e = row["jpl_e"]
        jpl_i = row["jpl_i"]
        jpl_Omega = row["jpl_loan"]
        jpl_omega = row["jpl_aop"]
        jpl_M = row["jpl_ma"]

        print("a: ", a, "       | JPL: ", jpl_a, " | Difference: ", abs((a-jpl_a)/jpl_a)*100, "%")
        print("e: ", e, "      | JPL: ", jpl_e, " | Difference: ", abs((e-jpl_e)/jpl_e)*100, "%")
        print("i: ", i, "      | JPL: ", jpl_i, " | Difference: ", abs((i-jpl_i)/jpl_i)*100, "%")
        print("Omega: ", Omega, "   | JPL: ", jpl_Omega, " | Difference: ", abs((Omega-jpl_Omega)/jpl_Omega)*100, "%")
        print("omega: ", omega, "  | JPL: ", jpl_omega, " | Difference: ", abs((omega-jpl_omega)/jpl_omega)*100, "%")
        print("M: ", M, "      | JPL: ", jpl_M, " | Difference: ", abs((M-jpl_M)/jpl_M)*100, "%")
        print("---------------------------------------------------------------------------------------------")

        

        
    
