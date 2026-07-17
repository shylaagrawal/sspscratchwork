# Shyla Agrawal
# 6/27
# computing "about equals"

def is_close_absolute(expected, actual, abs_tol):
    # Returns True if the absolute difference between "expected" and "actual" is within "abs_tol".
    # Otherwise, returns False.

    if (abs(expected - actual) < abs_tol):
        return True
    return False

def is_close_relative(expected, actual, rel_tol):
    # Returns True if the relative difference between "expected" and "actual" is within "rel_tol".
    # Otherwise, returns False
    
    if (abs(actual - expected) < abs(expected * rel_tol)):
        return True
    return False

def is_close(expected, actual, abs_tol, rel_tol):
    # Returns True if "expected" and "actual" are within an absolute difference "abs_tol",
    # OR within a relative difference "rel_tol". Otherwise, returns False.

    if (abs(expected - actual) < abs_tol):
        return True
    elif (abs(actual - expected) < abs(expected * rel_tol)):
        return True
    return False

def main():
    print('TESTING: is_close_absolute')
    assert is_close_absolute(1, 1, 0.1)
    assert not is_close_absolute(1, 2, 0.1)
    assert is_close_absolute(1, 2, 10)
    assert is_close_absolute(2, 1, 10)
    assert is_close_absolute(-100, -200, 1000)
    assert not is_close_absolute(-100, -200, 0.1)
    assert is_close_absolute(1.0, 1.0001, 0.01)
    print('All tests passed!')

    print('TESTING: is_close_relative')
    assert is_close_relative(1, 1.005, 0.01)
    assert not is_close_relative(1, 1.015, 0.01)
    assert is_close_relative(100, 101, 0.03)
    assert is_close_relative(100, 99, 0.03)
    assert not is_close_relative(100, 104, 0.03)
    assert not is_close_relative(100, 96, 0.03)
    assert not is_close_relative(0, 1, 0.1)
    print('All tests passed!')

    print('TESTING: is_close')
    # Within the absolute tolerance:
    assert is_close(10, 10.5, 2, 0.001)
    assert is_close(10, 9.5, 2, 0.001)
    assert is_close(0, -1.9999, 2, 0.001)
    assert is_close(-1.9999, 0, 2, 0.001)
    # Within the relative tolerance, but NOT the absolute tolerance:
    assert is_close(10, 10.5, 0.001, 0.1)
    assert is_close(10, 9.5, 0.001, 0.1)
    assert is_close(-10, -10.5, 0.001, 0.1)
    assert is_close(-10, -9.5, 0.001, 0.1)
    assert is_close(0.1, 0.105, 0.001, 0.1)
    # Within neither tolerance:
    assert not is_close(10, 12, 0.1, 0.1)
    assert not is_close(10, 8, 0.1, 0.1)
    assert not is_close(-10, -12, 0.1, 0.1)
    assert not is_close(-10, -8, 0.1, 0.1)
    assert not is_close(0, 1, 0.1, 0.1)
    print('All tests passed!')