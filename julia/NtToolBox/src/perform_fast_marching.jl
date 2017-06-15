using PyCall
@pyimport skfmm as fmm

function perform_fast_marching(D, start_points)
    """
        Implementation of the fast marching algorithm in 2Ds using the skfmm library
        D : 2D weight matrix, must be positive
        start_points : 2 x n array, start_points[:,i] is the ith starting point
    """
    D_temp = D
    D_temp[start_points[1,:] + (start_points[2,:] - 1).*size(D)[1]] = 0
    return fmm.distance(D_temp) + 1e-15
end
