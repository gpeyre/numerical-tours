import numpy as np
import skfmm as fmm

def perform_fast_marching(D,start_points):
    """
        Implementation of the fast marching algorithm in 2Ds using the skfmm library
        D : 2D weight matrix, must be positive
        start_points : 2 x n array, start_points[:,i] is the ith starting point
    """
    D_temp = np.copy(D)
    D_temp[start_points[0,:],start_points[1,:]] = 0
    return fmm.distance(D_temp) + 1e-15