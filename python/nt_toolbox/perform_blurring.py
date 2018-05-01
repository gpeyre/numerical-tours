import numpy as np
import pylab as pyl
from scipy import signal

def perform_blurring(M, sigma, bound="sym"):
    """
        perform_blurring - gaussian blurs an image

        M = perform_blurring(M, sigma, options);

        M is the original data
        sigma is the width of blurs (in pixels)

        Copyright (c) 2007 Gabriel Peyre
    """

    if np.all(sigma == 0):
        return M

    if np.ndim(M) > 2:
        for i in range(np.shape(M)[2]):
             M[:,:,i] = perform_blurring(M[:,:,i], sigma, bound);

    n = max(np.shape(M))

    eta = 4
    p = np.round((sigma*eta)/2.)*2+1
    p = np.minimum(p,(round(n/2.)*2-1)*np.ones(len(p)))

    A = np.array([1.,1.])
    if np.ndim(M) == 1:
        A = 1 #1D

    h = compute_gaussian_filter(p*A,sigma/(4.*n),n*A)
    M = perform_convolution(M, h, bound)

    return M


def compute_gaussian_filter(n,s,N):
    """
        compute_gaussian_filter - compute a 1D or 2D Gaussian filter.

          f = compute_gaussian_filter(n,s,N);

          'n' is the size of the filter, odd for no phase in the filter.
              (if too small it will alterate the filter).
              use n=[n1,n2] for a 2D filter or n = [n1] for a 1D filter
          's' is the standard deviation of the filter.
          'N' is the size of the big signal/image (supposed to lie in [0,1] or [0,1]x[0,1]).
              use N=[N1,N2] for a 2D filter or N = [N1] for a 1D filter

          The equation (in 1D) is
              f[k] = exp( -(x(k)^2/(2*s^2)) );
          where x spans [-1/2,1/2].

          The filter is normalised so that it sums to 1.

          Copyright (c) 2004 Gabriel Peyre
    """
    nd = 1;
    if len(n) > 1 and n[1] > 1:
        nd = 2

    if nd == 2 and len(s) == 1:
        s = np.hstack((s,s))

    if nd == 2 and len(N) == 1:
        N = np.hstack((N,N))

    if nd == 1 :
        f = build_gaussian_filter_1d(n,s,N)
    else:
        f = build_gaussian_filter_2d(n,s,N)
    return f

def build_gaussian_filter_2d(n,s,N=[]):
    """
        build_gaussian_filter_2d - compute a 2D Gaussian filter.

        f = build_gaussian_filter_2d(n,s,N);

        'n' is the size of the filter, odd for no phase in the filter.
            (if too small it will alterate the filter).
        's' is the standard deviation of the filter.
        'N' is the size of the big image (supposed to lie in [0,1]x[0,1]).

        The filter is normalised so that it sums to 1.

        Copyright (c) 2004 Gabriel Peyre
    """

    n = np.asarray(n)
    s = np.asarray(s)
    N = np.asarray(N)

    if len(N) == 0:
        N = n

    if len(N) == 1 or N[0] == 1:
        N = np.hstack((N,N))

    if len(s) == 1 or s[0] == 1:
        s = np.hstack((s,s))

    if len(s[s <= 0]) > 0:
        f = np.zeros(n)
        f[np.round((n-1)/2).astype(int)] = 1
        return f

    x = (np.arange(0,n[0])-(n[0]-1)/2.)/(N[0]-1)
    y = (np.arange(0,n[1])-(n[1]-1)/2.)/(N[1]-1)
    [Y,X] = np.meshgrid(y,x)
    f = np.exp(-(X**2/(2*s[0]**2)) - (Y**2/(2*s[1]**2)))
    f = f/np.sum(f)
    return f


def build_gaussian_filter_1d(n,s,N=[]):
    """
        build_gaussian_filter_1d - compute a Gaussian filter.

        f = build_gaussian_filter_1d(n,s,N);

        Copyright (c) 2004 Gabriel Peyre
    """
    if len(N) == 0:
        N = n

    n = n[0]
    s = s[0]
    N = N[0]

    if s <= 0:
        f = np.zeros(n)
        f[np.round((n-1)/2)] = 1
        return f

    x = (np.arange(0,n)-(n-1)/2.)/(N-1)
    f = np.exp(-x**2/(2*s**2))
    f = f/np.sum(f)
    return f


def perform_convolution(x,h,bound="sym"):
    """
        perform_convolution - compute convolution with centered filter.

        y = perform_convolution(x,h,bound);

        The filter 'h' is centred at 0 for odd
        length of the filter, and at 1/2 otherwise.

        This works either for 1D or 2D convolution.
        For 2D the matrix have to be square.

        'bound' is either 'per' (periodic extension)
        or 'sym' (symmetric extension).

        Copyright (c) 2004 Gabriel Peyre
    """

    if bound not in ["sym", "per"]:
        raise Exception('bound should be sym or per')

    if np.ndim(x) == 3 and np.shape(x)[2] < 4:
        #for color images
        y = x;
        for i in range(np.shape(x)[2]):
            y[:,:,i] = perform_convolution(x[:,:,i],h, bound)
        return y

    if np.ndim(x) == 3 and np.shape(x)[2] >= 4:
        raise Exception('Not yet implemented for 3D array, use smooth3 instead.')

    n = np.shape(x)
    p = np.shape(h)

    nd = np.ndim(x)

    if nd == 1:
        n = len(x)
        p = len(h)

    if bound == 'sym':

                #################################
        # symmetric boundary conditions #
        d1 = np.asarray(p).astype(int)/2  # padding before
        d2 = p - d1 - 1    			    # padding after

        if nd == 1:
        ################################# 1D #################################
            nx = len(x)
            xx = np.vstack((x[d1:-1:-1],x,x[nx-1:nx-d2-1:-1]))
            y = signal.convolve(xx,h)
            y = y[p:nx-p-1]

        elif nd == 2:
        ################################# 2D #################################
            #double symmetry
            nx,ny=np.shape(x)
            xx = x
            xx = np.vstack((xx[d1[0]:-1:-1,:], xx, xx[nx-1:nx-d2[0]-1:-1,:]))
            xx = np.hstack((xx[:,d1[1]:-1:-1], xx, xx[:,ny-1:ny-d2[1]-1:-1]))
            y = signal.convolve2d(xx,h,mode="same")
            y = y[(2*d1[0]):(2*d1[0]+n[0]+1), (2*d1[1]):(2*d1[1]+n[1]+1)]

    else:

        ################################
        # periodic boundary conditions #

        if p > n:
            raise Exception('h filter should be shorter than x.')
        n = np.asarray(n)
        p = np.asarray(p)
        d = np.floor((p-1)/2.)
        print(n)
        if nd == 1:
            h = np.vstack((h[d:],np.vstack((np.zeros(n-p),h[:d]))))
            y = np.real(pyl.ifft(pyl.fft(x)*pyl.fft(h)))
        else:
            h = np.vstack((h[d[0]:,:],np.vstack((np.zeros([n[0]-p[0],p[1]]),h[:(d[0]),:]))))
            h = np.hstack((h[:,d[1]:],np.hstack((np.zeros([n[0],n[1]-p[1]]),h[:,:(d[1])]))))
            y = np.real(pyl.ifft2(pyl.fft2(x)*pyl.fft2(h)))
    return y
