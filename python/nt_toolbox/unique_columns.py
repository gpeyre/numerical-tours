import numpy as np

def unique_columns(A):
    sorted_idx = np.lexsort((A[1,:],A[0,:]))
    A = A[:,sorted_idx]
    B = np.hstack((np.ones([2,1]),np.diff(A)))
    R = np.ravel(np.maximum(B[0,:],B[1:,]))
    return A[:,R != 0]
