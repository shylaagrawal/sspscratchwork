# Shyla Agrawal
# 7/17/2026
# SSP Orbit Determination Code 4

import numpy as np
import odlib
import pandas as pd
from importlib import reload

reload(odlib)

def main():
    data = pd.read_csv("f_and_g_test_cases.csv")

    for _, row in data.iterrows():
        r_vec = np.array([row["r2_x"], row["r2_y"], row["r2_z"]])
        v_vec = np.array([row["v2_x"], row["v2_y"], row["v2_z"]])
        tau = row["tau"]
        order_four = row["order4"]

        f, g = odlib.f_and_g(r_vec, v_vec, tau, order_four)

        print(f"Computed f: {f:.6f} Computed g: {g:.6f}")
        print(f"Expected f: {row["f_output"]} Expected g: {row["g_output"]}")
        print()
