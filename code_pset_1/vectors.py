# Shyla Agrawal
# 6/27
# doing vector math
import numpy as np

def vector_magnitude(vec):
    # Returns the length of the vector "vec" (list)
    sum_squares = 0
    for i in range(len(vec)):
        sum_squares += vec[i] ** 2
    return np.sqrt(sum_squares)

def dot_product(vec1, vec2):
    # Returns the dot product of "vec1" (list) and "vec2" (list)
    # These vectors must have the same number of elements

    # Prevents this function from running if vec1 and vec2 have a different number of elements
    assert len(vec1)==len(vec2), 'vec1 and vec2 must have the same length'

    result = 0
    for i, j in zip(vec1, vec2):
        result += i * j
    return result

def cross_product(vec1, vec2):
    # Returns the cross product of "vec1" (list) and "vec2" (list)
    # Each vector must have exactly 3 elements

    # Prevents this function from running if vec1 or vec2 have the wrong number of elements
    assert len(vec1)==3, 'vec1 must have exactly 3 elements'
    assert len(vec2)==3, 'vec2 must have exactly 3 elements'

    return [vec1[1]*vec2[2] - vec1[2]*vec2[1], vec1[2]*vec2[0] - vec1[0]*vec2[2], vec1[0]*vec2[1] - vec1[1]*vec2[0]]

def main():

    print('TESTING: vector_magnitude')
    print()
    print('Actual:  ', vector_magnitude([1]))
    print('Expected:', 1)
    print()
    print('Actual:  ', vector_magnitude([1, 1, -1, -1]))
    print('Expected:', 2)
    print()
    print('Actual:  ', vector_magnitude([10, -10]))
    print('Expected:', 14.142136)
    print()

    ####################

    print('TESTING: dot_product')
    print()
    print('Actual:  ', dot_product([], []))
    print('Expected:', 0)
    print()
    print('Actual:  ', dot_product([2,4], [1,1]))
    print('Expected:', 6)
    print()
    print('Actual:  ', dot_product([2,5,6], [3,7,8]))
    print('Expected:', 89)
    print()
    print('Actual:  ', dot_product([0,2], [2,0]))
    print('Expected:', 0)
    print()
    print('Actual:  ', dot_product([3], [6]))
    print('Expected:', 18)
    print()
    
    ####################
    
    print('TESTING: cross_product')
    print()
    print('Actual:  ', cross_product([1,0,0], [0,1,0]))
    print('Expected:', [0,0,1])
    print()
    print('Actual:  ', cross_product([2,5,6], [3,7,8]))
    print('Expected:', [-2,2,-1])
    print()
    print('Actual:  ', cross_product([1,5,9], [0,8,5]))
    print('Expected:', [-47, -5, 8])
    print()
    print('Actual:  ', cross_product([0,2,0], [2,0,2]))
    print('Expected:', [4, 0, -4])
    print()
