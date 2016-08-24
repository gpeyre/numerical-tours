import numpy as np
from scipy import linalg
from scipy import interpolate

def perform_image_similitude(M,mapping_type, u,u1,v,v1, w,w1):
    """
        perform_image_similitude
        
          M1 = perform_image_similitude(M,mapping_type,u,u1,v,v1,w,w1);
        
        Compute the affine similitude that map u to u1
        and v to v1, and then resample the image M.
        p and p1 are assumed to be in [0,1]≤
        
          If mapping_type=='similitude', compute a true similitude
              T(x,y) = [a -b] * [x] + [c]
                       [b  a]   [y]   [d]
          Solve the equations T(u)=u1 and T(v)=v1.
        
          If mapping_type=='similitude', compute a true similitude
              T(x,y) = [a  b] * [x] + [e]
                       [c  d]   [y]   [f]
          Solve the equations T(u)=u1 and T(v)=v1 and T(w)=w1.
        
        Copyright (c) 2006 Gabriel PeyrÈ
    """
    if mapping_type == "affine":   
        # the matrix of the linear system
        A = np.array([
             [u[0], u[1], 0   , 0   , 1, 0],
             [0   , 0   , u[0], u[1], 0, 1],
             [v[0], v[1], 0   , 0   , 1, 0],
             [0   , 0   , v[0], v[1], 0, 1],
             [w[0], w[1], 0   , 0   , 1, 0],
             [0   , 0   , w[0], w[1], 0, 1]])
        # the right hand size
        rhs = np.hstack((u1, v1, w1))
        # solve
        z = linalg.solve(A,rhs)
        # the similitude
        Q = np.array([[z[0], z[1]],[z[2], z[3]]])
        # the translation
        t = np.array([z[4], z[5]])
 
    else:
        raise Error('Unknown mapping');

    ### perform resampling ###
    
    # original grid in the warped domain
    n = np.shape(M)[0]
    x = np.linspace(0,1,n)
    [X,Y] = np.meshgrid(x,x)
    # inverse warping P1=T^-1(P)
    P = np.vstack((np.ravel(X),np.ravel(Y))) # position of the sampling
    P[0,:] = P[0,:] - t[0] # substract translation
    P[1,:] = P[1,:] - t[1]
    P1 = np.dot((np.linalg.inv(Q)),P) # undo similitude
    # reshape the results
    X1 = P1[0,:]
    Y1 = P1[1,:]
    
    M1 = np.reshape(interpolate.RectBivariateSpline(x,x,M).ev(X1,Y1), (n,n))
    return M1