import numpy as np

def div(Px,Py, bound="sym", order=1):
    """
        div - divergence operator

        fd = div(Px,Py, options);
        fd = div(P, options);

          options.bound = 'per' or 'sym'
          options.order = 1 (backward differences)
                        = 2 (centered differences)

          Note that the -div and grad operator are adjoint
          of each other such that 
              <grad(f),g>=<f,-div(g)>

          See also: grad.

        Copyright (c) 2007 Gabriel Peyre
    """
    
    # retrieve number of dimensions
    nbdims = np.ndim(Px)
    
    if nbdims >= 3:
        if nbdims == 3:
            Py = P[:,:,1]
            Px = P[:,:,0]
            nbdims=2
        else:
            Pz = P[:,:,:,2]
            Py = P[:,:,:,1]
            Px = P[:,:,:,0] 
            nbdims = 3
            
    if bound == "sym":
        nx = np.shape(Px)[0]
        if order == 1:
            fx = Px - Px[np.hstack(([0],np.arange(0,nx-1))),:]         
            fx[0,:]   = Px[0,:]                        # boundary
            fx[nx-1,:] = -Px[nx-2,:] 
            
            if nbdims >= 2:
                ny = np.shape(Py)[1]
                fy = Py - Py[:,np.hstack(([0],np.arange(0,ny-1)))]         
                fy[:,0]   = Py[:,0]                    # boundary
                fy[:,ny-1] = -Py[:,ny-2]  
                
            if nbdims >= 3:
                nz = np.shape(Pz)[2]
                fz = Pz - Pz[:,:,np.hstack(([0],np.arange(0,nz-1)))]         
                fz[:,:,0]   = Pz[:,:,0]                # boundary
                fz[:,:,nz-1] = -Pz[:,:,nz-2]        
        else:
            fx = (Px[np.hstack((np.arange(1,nx),[nx-1])),:] - Px[np.hstack(([0],np.arange(0,nx-1))),:])/2.
            fx[0,:] = + Px[1,:]/2. + Px[0,:]           # boundary
            fx[1,:] = + Px[2,:]/2. - Px[0,:]
            fx[nx-1,:] = - Px[nx-1,:]-Px[nx-2,:]/2.
            fx[nx-2,:] = + Px[nx-1,:]-Px[nx-3,:]/2.

            if nbdims >= 2:
                ny = np.shape(Py)[1]
                fy = (Py[:,np.hstack((np.arange(1,ny),[ny-1]))] - Py[:,np.hstack(([0],np.arange(0,ny-1)))])/2.
                fy[:,0] = + Py[:,1]/2. + Py[:,0]       # boundary
                fy[:,1] = + Py[:,2]/2. - Py[:,0]
                fy[:,ny-1] = - Py[:,ny-1]-Py[:,ny-2]/2.
                fy[:,ny-2] = + Py[:,ny-1]-Py[:,ny-3]/2.
            
            if nbdims >= 3:
                nz = np.shape(Pz)[2]
                fz = (Pz[:,:,np.hstack((np.arange(1,nz),[nz-1]))] - Pz[:,:,np.hstack(([0],np.arange(0,nz-1)))])/2.
                fz[:,:,0] = + Pz[:,:,1]/2. + Pz[:,:,0] # boundary
                fz[:,:,1] = + Pz[:,:,2]/2. - Pz[:,:,0]
                fz[:,:,ny-1] = - Pz[:,:,nz-1]-Pz[:,:,nz-2]/2.
                fz[:,:,ny-2] = + Pz[:,:,nz-1]-Pz[:,:,nz-3]/2.
    else:
        if order == 1:
            nx = np.shape(Px)[0]
            fx = Px-Px[np.hstack(([nx-1],np.arange(0,nx-1))),:]
            
            if nbdims >= 2:
                ny = np.shape(Py)[1]
                fy = Py-Py[:,np.hstack(([ny-1],np.arange(0,ny-1)))]
            
            if nbdims>=3:
                nz = np.shape(Pz)[2]
                fz = Pz-Pz[:,:,np.hstack(([nz-1],np.arange(0,nz-1)))]
            
        else:
            nx = np.shape(Px)[0]
            fx = (Px[np.hstack((np.arange(1,nx),[0])),:]) - (Px[np.hstack(([nx-1],np.arange(0,nx-1))),:])
            
            if nbdims >= 2:
                ny = np.shape(Py)[1]
                fy = (Py[:,np.hstack((np.arange(1,ny),[0]))]) - (Py[:,np.hstack(([ny-1],np.arange(0,ny-1)))])
            
            if nbdims >= 3:
                nz = np.shape(Pz)[2]
                fz = (Pz[:,:,np.hstack((np.arange(1,nz),[0]))]) - (Pz[:,:,np.hstack(([nz-1],np.arange(0,nz-1)))])

    # gather result
    if nbdims == 3:
        fd = fx+fy+fz
        
    elif nbdims == 2:
        fd = fx+fy

    else:
        fd = fx
        
    return fd