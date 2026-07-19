# Shyla Agrawal
# 7/18/2026
# SSP Orbit Determination Code 5

import numpy as np
import odlib
import pandas as pd
from importlib import reload

reload(odlib)

def main():
    data = pd.read_csv("mog_test_cases.csv")