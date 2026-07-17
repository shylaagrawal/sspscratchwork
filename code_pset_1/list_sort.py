# Shyla Agrawal
# 6/27
# manipulating lists and sorting them

def get_minimum(ls):
    smallest = ls[0]
    for i in ls:
        if i < smallest:
            smallest = i
    return smallest

def remove_minimum(ls):
    ls.remove(get_minimum(ls))
    return ls

def sort_list(ls):
    sorted = []
    for _ in range(len(ls)):
        sorted.append(get_minimum(ls))
        ls = remove_minimum(ls)
    return sorted

def main():

    print('TESTING: get_minimum')
    print()
    print('Actual:  ', get_minimum([1, 2, 3, 4, 5]))
    print('Expected:', 1)
    print()
    print('Actual:  ', get_minimum([6, 7, -1, 1]))
    print('Expected:', -1)
    print()
    print('Actual:  ', get_minimum([0, -10, 273849, 34]))
    print('Expected:', -10)
    print()
    print('Actual:  ', get_minimum([0, -3.5, -324, 9]))
    print('Expected:', -324)
    print()

    ####################

    print('TESTING: remove_minimum')
    print()
    print('Actual:  ', remove_minimum([2, 3, 4, 5, 6, 6, 7]))
    print('Expected:', [3, 4, 5, 6, 6, 7])
    print()
    print('Actual:  ', remove_minimum([2, 4]))
    print('Expected:,' [4])
    print()
    print('Actual:  ', remove_minimum([2, 0, -3, 5, 6]))
    print('Expected:', [2, 0, 5, 6])
    print()
    print('Actual:  ', remove_minimum([0, 3, 4, 5]))
    print('Expected:', [3, 4, 5])
    print()
    
    ####################
    
    print('TESTING: sort_list')
    print()
    print('Actual:  ', sort_list([0,1,0]))
    print('Expected:', [0,0,1])
    print()
    print('Actual:  ', sort_list([3,7,8]))
    print('Expected:', [3,7,8])
    print()
    print('Actual:  ', sort_list([0,8,5]))
    print('Expected:', [0,5,8])
    print()
    print('Actual:  ', sort_list([-2,0,2]))
    print('Expected:', [-2,0,2])
    print()
    