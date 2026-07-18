# Shyla Agrawal
# 7/17/2026
# SSP Orbit Determination Code 4

import numpy as np
import odlib
from importlib import reload
reload(odlib)

def main():
    data = np.loadtxt("agrawalOD4Input.txt")
    