# Shyla Agrawal
# 7/16/2026
# SSP Orbit Determination Code 3

import numpy as np
import odlib
from importlib import reload
reload(odlib)

def main():
    with open("agrawalOD3Input.txt") as f:
        lines = f.readlines()
    
    data = []
    for line in lines:
        line = line.strip()
        if line != "" and not line.startswith("#"):
            data.append(line)

    dec_positive = data[11].startswith("+")
    data[11] = data[11].replace("+", "").replace("-", "")

    data = [float(x) for x in data]

    jd_ref = data[0]
    orbital_elements = np.array(data[1:7])
    jd_target = data[7]
    sun_vec = np.array([data[14], data[15], data[16]])

    RA_pred, Dec_pred = odlib.ephemeris(*orbital_elements, jd_ref, jd_target, sun_vec)
    target_RA_degrees = odlib.hms_to_deg(data[8], data[9], data[10])
    target_dec_degrees = odlib.dms_to_deg(dec_positive, data[11], data[12], data[13])

    print("Predicted RA:", RA_pred)
    print("Expected RA:", target_RA_degrees)
    print()
    print("Predicted Dec:", Dec_pred)
    print("Expected Dec:", target_dec_degrees)


