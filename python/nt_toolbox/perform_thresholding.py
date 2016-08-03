import numpy as np

def perform_thresholding(f,M,type):
    """
        Only 3 types of thresholding currently implemented
    """
    if type == "largest":
        a = np.sort(np.ravel(abs(f)))[::-1] #sort a 1D copy of F in descending order
        T = a[M]
        y = f*(abs(f) > T)
    elif type == "soft":
        s = abs(f) - M
        s = (s + abs(s))/2
        y = np.sign(f)*s
    elif type == "hard":
        y = f*(abs(f) > M)
    return y