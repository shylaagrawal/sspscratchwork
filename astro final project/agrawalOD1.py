# Shyla Agrawal
# 7/8/2026
# SSP Orbit Determination Code 1

import numpy as np
import odlib
from importlib import reload
reload(odlib)

def main():
    data = np.loadtxt("agrawalOD1Input.txt")
    r_vec = data[0]
    v_vec = data[1]
    h = odlib.angular_momentum(r_vec, v_vec)
    print(f"{h[0]:.6f}, {h[1]:.6f}, {h[2]:.6f}")