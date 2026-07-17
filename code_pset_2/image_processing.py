# Shyla Agrawal
# 6/29
# library that manipulates images in various ways

import matplotlib.pyplot as plt
import numpy as np

def create_inverted(image_filename):
    # Reads the image ("image_filename") and returns a color-inverted version

    # Reads the image
    image = plt.imread(image_filename)
    return 1 - image[:,:,:]
    

def create_grayscale(image_filename):
    # Reads the image ("image_filename") and returns it in grayscale
        
    # Reads the image
    image = plt.imread(image_filename)
    for row in range(len(image)):
        for col in range(len(image[row])):
            image[row, col] = (image[row,col,0] + image[row,col,1] + image[row,col,2])/3

    return image

def create_bluescreen(image_filename, background_filename):
    # Reads the main image ("image_filename") and the background image ("background_filename"),
    # and returns the main image with every blue pixel replaced by the background image
    image = plt.imread(image_filename)
    background = plt.imread(background_filename)
    for row in range(len(image)):
        for col in range(len(image[row])):
            pixel_colors = image[row, col]
            if pixel_colors[2] > pixel_colors[1] + pixel_colors[0]:
                image[row, col] = background[row, col]
               
    return image
