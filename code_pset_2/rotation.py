# Shyla Agrawl
# 6/29
# library that rotates a vector based on an angle

import numpy as np

def rotate_2d(vec2d, theta_rad):
    # Returns the vector "vec2d" rotated by "theta_rad" radians

    cos = np.cos(theta_rad)
    sin = np.sin(theta_rad)
    rotation_matrix = [[cos, -sin], 
                       [sin, cos]]

    return np.round(rotation_matrix @ vec2d, decimals = 8)


def print_test_result(actual, expected):
    # A helper function for printing out tests. Prints "actual", "expected", and "actual-expected"
    print('Actual:  ', actual)
    print('Expected:', expected)
    print('Difference:', actual-expected)
    if np.any(np.abs(actual-expected) > 1e-4):
        print('ERROR: Difference is too large!')
    print()

def main():
    
    print('TESTING: rotate_2d')
    print()

    # Rotates the x-unit vector [1, 0]
    print_test_result(rotate_2d(np.array([1,0]), 0), np.array([1,0]))
    print_test_result(rotate_2d(np.array([1,0]), np.pi/2), np.array([0,1]))
    print_test_result(rotate_2d(np.array([1,0]), np.pi), np.array([-1,0]))
    
    # Rotates the y-unit vector [0, 1]
    print_test_result(rotate_2d(np.array([0,1]), 0), np.array([0,1]))
    print_test_result(rotate_2d(np.array([0,1]), np.pi/2), np.array([-1,0]))

    # Rotates more complicated vectors
    print_test_result(
        rotate_2d(np.array([3,7]), 1.7*np.pi),
        np.array([7.42647472, 1.68744578])
    )
    print_test_result(
        rotate_2d(np.array([-100, 300]), -0.1*np.pi),
        np.array([-2.40055332, 316.21865433])
    )

