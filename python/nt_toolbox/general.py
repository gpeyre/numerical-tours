import numpy as np
import pylab
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from scipy import ndimage

# general.py TODO: try to not make use of transform.resize
from skimage import transform  ## commented


def crop(M, n=None):
    """Return a centered square crop, including the full-size case."""
    M = np.asarray(M)
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("crop expects a square grayscale image")
    n = M.shape[0] // 2 if n is None else int(n)
    if not 1 <= n <= M.shape[0]:
        raise ValueError("crop size must be between 1 and the image size")
    start = (M.shape[0] - n) // 2
    return M[start : start + n, start : start + n]


def circshift(x, p):
    """Circularly shift toward lower indices, matching circshift1d."""
    p = tuple(int(k) for k in p)
    return np.roll(x, tuple(-k for k in p), axis=tuple(range(len(p))))


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


def rescale(f, a=0, b=1):
    """
    Rescale linearly the dynamic of a vector to fit within a range [a,b]
    """
    v = f.max() - f.min()
    g = (f - f.min()).copy()
    if v > 0:
        g = g / v
    return a + g * (b - a)


def reverse(x):
    """
    Reverse a vector.
    """
    return x[::-1]
