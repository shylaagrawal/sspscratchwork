# Shyla Agrawal
# 7/16/2026
# SSP Orbit Determination Code 2

import numpy as np
import odlib
from importlib import reload
reload(odlib)

def main():
    data = np.loadtxt("agrawalOD2Input.txt")
    jd_ref = data[0]

    a = data[1]
    e = data[2]
    i = data[3]
    Omega = data[4]
    omega = data[5]
    M = data[6]

    jd_target = data[7]

    RA_h = data[8]
    RA_m = data[9]
    RA_s = data[10]

    Dec_d = data[11]
    Dec_m = data[12]
    Dec_s = data[13]
    
    sun_vec = np.array([data[14], data[15], data[16]])

    