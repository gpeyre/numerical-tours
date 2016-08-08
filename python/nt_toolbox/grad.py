import numpy as np

def grad(M, bound="sym", order=1):
    """
        grad - gradient, forward differences
        
          [gx,gy] = grad(M, options);
        or
          g = grad(M, options);
        
          options.bound = 'per' or 'sym'
          options.order = 1 (backward differences)
                        = 2 (centered differences)
        
          Works also for 3D array.
          Assme that the function is evenly sampled with sampling step 1.
        
          See also: div.
        
          Copyright (c) Gabriel Peyre
    """    


    # retrieve number of dimensions
    nbdims = np.ndim(M)
    
    
    if bound == "sym":  
        nx = np.shape(M)[0]
        if order == 1:
            fx = M[np.hstack((np.arange(1,nx),[nx-1])),:] - M
        else:
            fx = (M[np.hstack((np.arange(1,nx),[nx-1])),:] - M[np.hstack(([0],np.arange(0,nx-1))),:])/2.
            # boundary
            fx[0,:] = M[1,:]-M[0,:]
            fx[nx-1,:] = M[nx-1,:]-M[nx-2,:]
            
        if nbdims >= 2:
            ny = np.shape(M)[1]
            if order == 1:
                fy = M[:,np.hstack((np.arange(1,ny),[ny-1]))] - M
            else:
                fy = (M[:,np.hstack((np.arange(1,ny),[ny-1]))] - M[:,np.hstack(([0],np.arange(ny-1)))])/2.
                # boundary
                fy[:,0] = M[:,1]-M[:,0]
                fy[:,ny-1] = M[:,ny-1]-M[:,ny-2]
    
        if nbdims >= 3:
            nz = np.shape(M)[2]
            if order == 1:
                fz = M[:,:,np.hstack((np.arange(1,nz),[nz-1]))] - M
            else:
                fz = (M[:,:,np.hstack((np.arange(1,nz),[nz-1]))] - M[:,:,np.hstack(([0],np.arange(nz-1)))])/2.
                # boundary
                fz[:,:,0] = M[:,:,1]-M[:,:,0]
                fz[:,:,ny-1] = M[:,:,nz-1]-M[:,:,nz-2]            
    else:
        nx = np.shape(M)[0]
        if order == 1:
            fx = M[np.hstack((np.arange(1,nx),[0])),:] - M
        else:
            fx = (M[np.hstack((np.arange(1,nx),[0])),:] - M[np.hstack(([nx-1],np.arange(nx-1))),:])/2.
            
        if nbdims >= 2:
            ny = np.shape(M)[1]
            if order == 1:
                fy = M[:,np.hstack((np.arange(1,ny),[0]))] - M
            else:
                fy = (M[:,np.hstack((np.arange(1,ny),[0]))] - M[:,np.hstack(([ny-1],np.arange(ny-1)))])/2.
        
        if nbdims >= 3:
            nz = np.shape(M)[2]
            if order == 1:
                fz = M[:,:,np.hstack((np.arange(1,nz),[0]))] - M
            else:
                fz = (M[:,:,np.hstack((np.arange(1,nz),[0]))] - M[:,:,np.hstack(([nz-1],np.arange(nz-1)))])/2.   
   
    if nbdims==2:
        fx = np.concatenate((fx[:,:,np.newaxis],fy[:,:,np.newaxis]), axis=2)
    elif nbdims==3:
        fx = np.concatenate((fx[:,:,:,np.newaxis],fy[:,:,:,np.newaxis],fz[:,:,:,np.newaxis]),axis=3)
    
    return fx