import numpy as np
import pylab
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from scipy import ndimage
# TODO: try to not make use of transform.resize
from skimage import transform

def circshift(x, p):
    """
        Circular shift of an array.
    """
    y = x.copy()
    y = np.concatenate((y[p[0]::, :], y[:p[0]:, :]), axis=0)
    if x.shape[1] > 0 and len(p) > 1:
        y = np.concatenate((y[:, p[0]::], y[:, :p[0]:]), axis=1)
    return y

def circshift1d(x, k):
    """ 
        Circularly shift a 1D vector
    """
    return np.roll(x, -k, axis=0)

def clamp(x, a=[], b=[]):
    """
     clamp - clamp a value

       y = clamp(x,a,b);

     Default is [a,b]=[0,1].

       Copyright (c) 2004 Gabriel Peyre
    """

    if a == []:
        a = 0.0
    if b == []:
        b = 1.0
    return np.minimum(np.maximum(x, a), b)

def rescale(f,a=0,b=1):
    """
        Rescale linearly the dynamic of a vector to fit within a range [a,b]
    """
    v = f.max() - f.min()
    g = (f - f.min()).copy()
    if v > 0:
        g = g / v
    return a + g*(b-a)

def reverse(x):
    """
        Reverse a vector. 
    """
    return x[::-1]