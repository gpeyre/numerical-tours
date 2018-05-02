import numpy as np

def perform_wavelet_transf(f, Jmin, dir, filter = "9-7",separable = 0, ti = 0):

    """""
    perform_wavelet_transf - peform fast lifting transform

    y = perform_wavelet_transf(x, Jmin, dir, filter = "9-7",separable = 0, ti = 0);

    Implement 1D and 2D symmetric wavelets with symmetric boundary treatements, using
    a lifting implementation.

    filter gives the coefficients of the lifting filter.
    You can use h='linear' or h='7-9' to select automatically biorthogonal
    transform with 2 and 4 vanishing moments.

    You can set ti=1 to compute a translation invariant wavelet transform.

    You can set separable=1 to compute a separable 2D wavelet
    transform.

    Copyright (c) 2008 Gabriel Peyre
    """

    #copy f
    x = np.copy(f)

    #convert Jmin to int
    Jmin = int(Jmin)

    # detect dimensionality
    d = np.ndim(x)
    # P/U/P/U/etc the last coefficient is scaling
    if filter in ["linear","5-3"]:
        h = [1/2, 1/4, np.sqrt(2)]

    elif filter in ["9-7","7-9"]:
        h = [1.586134342, -.05298011854, -.8829110762, .4435068522, 1.149604398]

    else:
        raise ValueError('Unknown filter')

    if d == 2 and separable == 1:
        ti = 0
        if ti == 1:
            wrn.warning("Separable does not works for translation invariant transform")

        # perform a separable wavelet transform
        n = np.shape(x)[0]
        if dir == 1:
            for i in range(n):
                x[:,i] = perform_wavelet_transf(x[:,i], Jmin, dir, filter, separable, ti)
            for i in range(n):
                x[i,:] = np.transpose(perform_wavelet_transf(np.transpose(x[i,:]), Jmin, dir, filter, separable, ti))
        else:
            for i in range(n):
                x[i,:] = np.transpose(perform_wavelet_transf(np.transpose(x[i,:]), Jmin, dir, filter, separable, ti))
            for i in range(n):
                x[:,i] = perform_wavelet_transf(x[:,i], Jmin, dir, filter, separable, ti)


    # number of lifting steps
    if np.ndim(x) == 1:
        n = len(x)
    else:
        n = np.shape(x)[1]
    m = (len(h)-1)//2
    Jmax = int(np.log2(n)-1)
    jlist = range(Jmax,Jmin-1,-1)

    if dir == -1:
        jlist = range(Jmin,Jmax+1,1)

    if ti == 0:
        # subsampled
        for j in jlist:
            if d == 1:
                x[:2**(j+1),:] = lifting_step(x[:2**(j+1)], h, dir)
            else:
                x[:2**(j+1),:2**(j+1)] = lifting_step(x[:2**(j+1),:2**(j+1)], h, dir)
                x[:2**(j+1),:2**(j+1)] = np.transpose(lifting_step(np.transpose(x[:2**(j+1),:2**(j+1)]), h, dir))

    else:
        # TI
        nJ = Jmax - Jmin + 1
        if dir == 1 and d == 1:
            x = np.tile(x,(nJ + 1,1,1))
        elif dir == 1 and d == 2:
            x = np.tile(x,(3*nJ + 1,1,1))
        #elif dir == 1:
        #    x = np.tile(x,(1,1,1))
        for j in jlist:
            dist = 2**(Jmax - j)

            if d == 1:
                if dir == 1:
                    x[:(j-Jmin+2),:,:] = lifting_step_ti(x[0,:,:], h, dir, dist)
                else:
                    x[0,:,:] = lifting_step_ti(x[:(j-Jmin+2),:,:], h, dir, dist)
            else:
                dj = 3*(j-Jmin)

                if dir == 1:
                    x[[0,dj+1],:,:] = lifting_step_ti(x[0,:,:], h, dir, dist)

                    x[[0,dj+2],:,:] = lifting_step_ti(np.transpose(x[0,:,:]), h, dir, dist)
                    x[0,:,:] = np.transpose(x[0,:,:])
                    x[dj+2,:,:] = np.transpose(x[dj+2,:,:])

                    x[[1+dj,3+dj],:,:] = lifting_step_ti(np.transpose(x[dj+1,:,:]), h, dir, dist)
                    x[dj+1,:,:] = np.transpose(x[dj+1,:,:])
                    x[dj+3,:,:] = np.transpose(x[dj+3,:,:])
                else:

                    x[dj+1,:,:] = np.transpose(x[dj+1,:,:])
                    x[dj+3,:,:] = np.transpose(x[dj+3,:,:])

                    x[dj+1,:,:] = np.transpose(lifting_step_ti(x[[1+dj, 3+dj],:,:], h, dir, dist))

                    x[0,:,:] = np.transpose(x[0,:,:])
                    x[dj+2,:,:] = np.transpose(x[dj+2,:,:])
                    x[0,:,:] = np.transpose(lifting_step_ti(x[[0,dj+2],:,:], h, dir, dist))

                    x[0,:,:] = lifting_step_ti(x[[0,dj+1],:,:], h, dir, dist)

        if dir == -1:
            x = x[0,:,:]

    return x

###########################################################################
###########################################################################
###########################################################################

def lifting_step(x0, h, dir):

    #copy x
    x = np.copy(x0)

    # number of lifting steps
    m = (len(h) - 1)//2

    if dir==1:
        # split
        d = x[1::2,]
        x = x[0::2,]
        for i in range(m):
            d = d - h[2*i] * (x + np.vstack((x[1:,],x[-1,])))
            x = x + h[2*i+1] * (d + np.vstack((d[0,],d[:-1,])))
        x = np.vstack((x*h[-1],d/h[-1]))

    else:
        # retrieve detail coefs
        end = len(x)
        d = x[end//2:,]*h[-1]
        x = x[:end//2,]/h[-1]
        for i in range(m,0,-1):
            x = x - h[2*i-1] * (d + np.vstack((d[0,],d[:-1,])))
            d = d + h[2*i-2] * (x + np.vstack((x[1:,],x[-1,])))
        # merge
        x1 = np.vstack((x,x))
        x1[::2,] = x
        x1[1::2,] = d
        x = x1

    return x

###########################################################################
###########################################################################
###########################################################################
def lifting_step_ti(x0, h, dir, dist):

    #copy x
    x = np.copy(x0)

    # number of lifting steps
    m = (len(h) - 1)//2
    n = np.shape(x[0])[0]

    s1 = np.arange(1, n+1) + dist
    s2 = np.arange(1, n+1) - dist

    # boundary conditions
    s1[s1 > n] = 2*n - s1[s1 > n]
    s1[s1 < 1] = 2   - s1[s1 < 1]

    s2[s2 > n] = 2*n - s2[s2 > n]
    s2[s2 < 1] = 2   - s2[s2 < 1]

    #indices in python start from 0
    s1 = s1 - 1
    s2 = s2 - 1

    if dir == 1:
        # split
        d = x
        for i in range(m):
            if np.ndim(x) == 2 :
                x = np.tile(x,(1,1,1))
            d = d - h[2*i]   * (x[:,s1,:] + x[:,s2,:])
            x = x + h[2*i+1] * (d[:,s1,:] + d[:,s2,:])

        #merge
        x = np.concatenate((x*h[-1],d/h[-1]))

    else:
        # retrieve detail coefs

        d = x[1,:,:]*h[-1]
        x = x[0,:,:]/h[-1]

        for i in range(m,0,-1):
            x = x - h[2*i-1] * (d[s1,:] + d[s2,:])
            d = d + h[2*i-2] * (x[s1,:] + x[s2,:])

        # merge
        x = (x + d)/2

    return x
