import numpy as np
from nt_toolbox.perform_fast_marching import *


def perform_redistancing(D):
    """
        perform_redistancing - redistance a function

          D1 = perform_redistancing(D, options);

          Compute a signed distance function D1 that has the same 0 level set as D (2D-matrix).

          Copyright (c) 2007 Gabriel Peyre
    """

    n,p = np.shape(D)
    eps = np.finfo(float).eps
    
    # horizontal
    P1 = D[0:(n-1),:]
    P2 = D[1:,:]
    P = ((P1*P2) <= 0)
    
    d = abs(P1-P2)
    d[d < eps] = 1
    v1 = abs(P1)/d
    v2 = abs(P2)/d
    Ah = ((np.vstack((P,np.zeros([1,n]))) + np.vstack((np.zeros([1,n]),P))) > 0)
    Vh = np.maximum(np.vstack((v1,np.zeros([1,n]))), np.vstack((np.zeros([1,n]),v2)))
    
    # vertical
    P1 = D[:,0:(p-1)]
    P2 = D[:,1:]
    P = ((P1*P2) <= 0)
    
    d = abs(P1-P2)
    d[d < eps] = 1
    v1 = abs(P1)/d
    v2 = abs(P2)/d
    Av = ((np.hstack((P,np.zeros([p,1]))) + np.hstack((np.zeros([p,1]),P))) > 0)
    Vv = np.maximum(np.hstack((v1,np.zeros([p,1]))), np.hstack((np.zeros([p,1]),v2)))
    
    V = np.zeros([n,p])
    I = np.where(Ah > 0)
    V[I] = Vh[I]
    I = np.where(Av > 0)
    V[I] = np.maximum(V[I],Vv[I])
    
    I = np.where(V != 0)
    x,y = I[0],I[1]
    start_points = np.vstack((x,y))
    D1 = perform_fast_marching(np.ones([n,p]), start_points)
    D1 = D1*n
    D1[D < 0] = -D1[D < 0]
    
    return D1