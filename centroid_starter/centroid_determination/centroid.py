# Shyla Agrawal
# 7/4/26
# Gets centriod of a star in an image

import numpy as np
import matplotlib.pyplot as plt

def get_centroid(grid):
    # Computes the intensity-weighted centroid of an image
    grid = np.array(grid, dtype=float)
    rows, cols = grid.shape
    total_mass = np.sum(grid)
    weighted_x = weighted_y = 0

    for row_index in range(rows):
        weighted_x += row_index * np.sum(grid[row_index, :])

    for col_index in range(cols):
        weighted_y += col_index * np.sum(grid[:, col_index])

    centroid_x = weighted_x / total_mass
    centroid_y = weighted_y / total_mass

    return (np.round(centroid_x, decimals = 3), np.round(centroid_y, decimals = 3))

def get_subgrid(image_grid, row_center, col_center, grid_size):
    # Returns a subgrid of "image_grid". Its center should be "row_center", "col_center", and its side length should be "grid_size"
    # "grid_size" must be positive & odd

    half_length = int(grid_size/2)
    return image_grid[row_center - half_length : row_center + half_length + 1, col_center - half_length : col_center + half_length + 1]
    

def get_centroid_from_image(image_grid, row_center, col_center, grid_size):
    # Given "image_grid", an array of pixels representing an image, computes the centroid at ("row_center", "col_center") using a subgrid of size "grid_size"
    # "grid_size" must be positive & odd

    centroid = get_centroid(get_subgrid(image_grid, row_center, col_center, grid_size))
    return (centroid[0] + row_center - grid_size // 2, centroid[1] + col_center - grid_size // 2)
    

def main():

    # Small grid for testing get_centroid
    # The centroid is approximately: row = 2.191 and column = 1.824 
    test_grid_small = np.array([
        [0, 33, 21, 33, 8],
        [0, 56, 51, 53, 26],
        [23, 120, 149, 73, 18],
        [55, 101, 116, 50, 16],
        [11, 78, 26, 2, 10],
    ])

    print('--- Testing: get_centroid ---')

    values = get_centroid(test_grid_small)
    print("expected: row =", 2.191, "col =", 1.824)
    print("actual:   row =", values[0], "col =", values[1])
    
    ########################################

    # Larger grid for testing get_subgrid and get_centroid_from_image
    # All values are 0 except for the one at row=30, col=70
    test_image = np.zeros([100, 100])
    test_image[30, 70] = 1
    
    print('--- Testing: get_subgrid ---')

    subgrid = get_subgrid(test_image, 30, 70, 5)
    plt.imshow(subgrid)
    plt.show()

    ########################################

    print('--- Testing: get_centroid_from_image ---')

    values = get_centroid_from_image(test_image, 29, 71, 5)

    print("Expected: row = 30, col = 70")
    print("Actual:   row =", values[0], "col =", values[1])    
