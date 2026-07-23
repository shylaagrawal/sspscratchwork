# Shyla Agrawal
# 7/9/26
# contains LSPR methods

import numpy as np

def xy_to_radec(x_pix, y_pix, b1, a11, a12, b2, a21, a22):
    # Converts a location in an image (x_pix, y_pix) to an RA and DEC coordinate in degrees
    # Uses the 6 plate solve coefficients to do so
    
    return (b1 + a11 * x_pix + a12 * y_pix, b2 + a21 * x_pix + a22 * y_pix)

def print_test_result(actual, expected):
    # A helper function for printing out tests. Prints "actual", "expected", and "actual-expected"
    # If any differences are larger than 1e-4, prints an error message
    print('Actual:  ', actual)
    print('Expected:', expected)
    print('Difference:', actual-expected)
    if np.any(np.abs(actual-expected) > 1e-4):
        print('ERROR: Difference is too large!')
    print()

def main():

    # Test plate solve coefficients
    b1, a11, a12, b2, a21, a22 = (271.23586782, -0.00201951, -0.002526924, 14.692310024, 0.002779617, -0.001835918)
    
    # -- TEST 1 --
    x, y = 0, 0
    actual_ra, actual_dec = xy_to_radec(x, y, b1, a11, a12, b2, a21, a22)
    expected_ra, expected_dec = b1, b2
    print('TEST 1: x, y =', x, y)
    print()
    print_test_result(actual_ra, expected_ra)
    print_test_result(actual_dec, expected_dec)

    # -- TEST 2 --
    x, y = 100, 0
    actual_ra, actual_dec = xy_to_radec(x, y, b1, a11, a12, b2, a21, a22)
    expected_ra, expected_dec = 271.03391682, 14.970271724
    print('TEST 2: x, y =', x, y)
    print()
    print_test_result(actual_ra, expected_ra)
    print_test_result(actual_dec, expected_dec)

    # -- TEST 3 --
    x, y = 57.3, 201.1
    actual_ra, actual_dec = xy_to_radec(x, y, b1, a11, a12, b2, a21, a22)
    expected_ra, expected_dec = 270.6119854806, 14.4823789683
    print('TEST 3: x, y =', x, y)
    print()
    print_test_result(actual_ra, expected_ra)
    print_test_result(actual_dec, expected_dec)
