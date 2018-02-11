import numpy as np
import pylab as pyl

def perform_stft(x,w,q,n):
    """
        perform_stft - compute a local Fourier transform
        
        Forward transform:
          MF = perform_stft(M,w,q, options);
        Backward transform:
          M  = perform_stft(MF,w,q, options);
        
          w is the width of the window used to perform local computation.
          q is the spacing betwen each window.
        
          MF(:,i) contains the spectrum around point (i-1)*q
        
          A typical use, for an redundancy of 2 could be w=2*q+1
        
          options.bound can be either 'per' or 'sym'
          
          No multichannel, no Gabor stft
        
          options.normalization can be set to
              'tightframe': tight frame transform, with energy conservation.
              'unit': unit norm basis vectors, usefull to do thresholding
        
          Copyright (c) 2008 Gabriel Peyre
    """

    #bound = getoptions(options, 'bound', 'per');
    #transform_type = getoptions(options, 'transform_type', 'fourier');
    #normalization = getoptions(options, 'normalization', 'tightframe');
    #window_type = getoptions(options, 'window_type', 'sin');
    #eta = getoptions(options, 'eta', 1);
    
    if np.ndim(x) == 1:
        dir = 1
    else:
        dir = -1
    
    # perform sampling
    X = np.arange(1,n+2,q)

    p = len(X)
    eta = 1
    
    if w%2 == 1:
        w = np.ceil((w-1)/2)*2+1
        w1 = (w-1)//2
        dX = np.arange(-w1,w1+1)
    else:
        dX = np.arange(-w//2+1,w//2+1)
    
    X1 = np.tile(X,(w,1)) + np.transpose(np.tile(dX, (p,1)))
    #periodic boundary conditions
    X1 = ((X1-1)%n)+1;
    
    I = X1 -1
    
    # build a sin weight function
    W = .5 *(1 - np.cos(2*np.pi*np.arange(0,w)/(w-1)))
    
    #renormalize the windows
    weight = np.zeros(n)
    
    for i in range(p):
        weight[I[:,i]] = weight[I[:,i]] + W**2
    
    weight = np.sqrt(weight)
    Weight = np.transpose(np.tile(W, (p,1)))
    
    for i in range(p):
        Weight[:,i] = Weight[:,i]/weight[I[:,i]]
        
    #compute the transform
    if dir == 1:
        y = np.zeros([eta*w,p])
        if w%2 == 1:
            m = (eta*w+1)//2
            w1 = (w-1)//2
            sel = np.arange(m-w1,m+w1+1) - 1
        else:
            m = (eta*w)//2+1 
            w1 = w//2
            sel = np.arange(m-w1,m+w1) - 1
        y[sel,:] = x[I]*Weight

        #perform the transform
        y = my_transform(y,+1)

    else:
        x = my_transform(x,-1)
        x = np.real(x*Weight)
        y = np.zeros(n)
        for i in range(p):
            y[I[:,i]] = y[I[:,i]] + x[:,i]

    return y
########################################################################
def my_transform(x,dir):

    # my_transform - perform either FFT with energy conservation.
    # Works on array of size (w,w,a,b) on the 2 first dimensions.
    w = np.shape(x)[0]
    if dir == 1 :
        y = np.transpose(pyl.fft(np.transpose(x)))/np.sqrt(w)
    else :
        y = np.transpose(pyl.ifft(np.transpose(x)*np.sqrt(w)))
    return y