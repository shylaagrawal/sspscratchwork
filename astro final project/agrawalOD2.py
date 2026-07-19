# Shyla Agrawal
# 7/13/2026
# SSP Orbit Determination Code 2

import numpy as np
import odlib
from importlib import reload
reload(odlib)

def main():
    data = np.loadtxt("agrawalOD2Input.txt")
    r_vec = np.array([data[1], data[2], data[3]])
    v_vec = np.array([data[4], data[5], data[6]])
    a, e, i, Omega, omega, M = odlib.orbital_elements(r_vec, v_vec)
    print(f"a:       {a:.6f}           expected: {data[7]:.6f}")
    print(f"e:       {e:.6f}           expected: {data[8]:.6f}")
    print(f"i:       {i:.6f}          expected: {data[9]:.6f}")
    print(f"Omega:   {Omega:.6f}         expected: {data[10]:.6f}")
    print(f"omega:   {omega:.6f}         expected: {data[11]:.6f}")
    print(f"M:       {M:.6f}         expected: {data[12]:.6f}")