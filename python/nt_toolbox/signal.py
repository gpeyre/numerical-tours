import numpy as np
import pylab
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from scipy import ndimage
# signal.py: TODO: try to not make use of transform.resize
from skimage import transform

from . import general as nt

#from nt_toolbox.general import *

def bilinear_interpolate(im, x, y):
    x = np.asarray(x)
    y = np.asarray(y)

    x0 = np.floor(x).astype(int)
    x1 = x0 + 1
    y0 = np.floor(y).astype(int)
    y1 = y0 + 1

    x0 = np.clip(x0, 0, im.shape[1]-1);
    x1 = np.clip(x1, 0, im.shape[1]-1);
    y0 = np.clip(y0, 0, im.shape[0]-1);
    y1 = np.clip(y1, 0, im.shape[0]-1);

    Ia = im[ y0, x0 ]
    Ib = im[ y1, x0 ]
    Ic = im[ y0, x1 ]
    Id = im[ y1, x1 ]

    wa = (x1-x) * (y1-y)
    wb = (x1-x) * (y-y0)
    wc = (x-x0) * (y1-y)
    wd = (x-x0) * (y-y0)

    return wa*Ia + wb*Ib + wc*Ic + wd*Id

def cconv(x, h, d):
    """
        Circular convolution along dimension d.
        h should be small and with odd size
    """
    if d == 2:
        # apply to transposed matrix
        return np.transpose(cconv(np.transpose(x), h, 1))
    y = np.zeros(x.shape)
    p = len(h)
    pc = int(round( float((p - 1) / 2 )))
    for i in range(0, p):
        y = y + h[i] * nt.circshift1d(x, i - pc)
    return y

def div(g):
    """
        Compute a finite difference approximation of the gradient of a 2D vector field, assuming periodic BC.
    """
    S = g.shape;
    s0 = np.concatenate( ([S[0]-1], np.arange(0,S[0]-1)) )
    s1 = np.concatenate( ([S[1]-1], np.arange(0,S[1]-1)) )
    f = (g[:,:,0] - g[s0,:,0]) + (g[:,:,1]-g[:,s1,1])
    return f

def gaussian_blur(f, sigma):

    """ gaussian_blur - gaussian blurs an image
    %
    %   M = perform_blurring(M, sigma, options);
    %
    %   M is the original data
    %   sigma is the std of the Gaussian blur (in pixels)
    %
    %   Copyright (c) 2007 Gabriel Peyre
    """
    if sigma<=0:
        return;
    n = max(f.shape);
    t = np.concatenate( (np.arange(0,n/2+1), np.arange(-n/2,-1)) )
    [Y,X] = np.meshgrid(t,t)
    h = np.exp( -(X**2+Y**2)/(2.0*float(sigma)**2) )
    h = h/np.sum(h)
    return np.real( pylab.ifft2(pylab.fft2(f) * pylab.fft2(h)) )

def grad(f):
    """
        Compute a finite difference approximation of the gradient of a 2D image, assuming periodic BC.
    """
    S = f.shape;
#   g = np.zeros([n[0], n[1], 2]);
    s0 = np.concatenate( (np.arange(1,S[0]),[0]) )
    s1 = np.concatenate( (np.arange(1,S[1]),[0]) )
    g = np.dstack( (f[s0,:] - f, f[:,s1] - f))
    return g

def imageplot(f, str='', sbpt=[]):
    """
        Use nearest neighbor interpolation for the display.
    """
    if sbpt != []:
        plt.subplot(sbpt[0], sbpt[1], sbpt[2])
    imgplot = plt.imshow(f, interpolation='nearest')
    imgplot.set_cmap('gray')
    pylab.axis('off')
    if str != '':
        plt.title(str)

def load_image(name, n=-1, flatten=1, resc=1, grayscale=1):
    """
        Load an image from a file, rescale its dynamic to [0,1], turn it into a grayscale image
        and resize it to size n x n.
    """
    f = plt.imread(name)
    # turn into normalized grayscale image
    if grayscale == 1:
        if (flatten==1) and np.ndim(f)>2:
            f = np.sum(f, axis=2)
    if resc==1:
        f = nt.rescale(f)
    # change the size of the image
    if n > 0:
        if np.ndim(f)==2:
            f = transform.resize(f, [n, n], 1)
        elif np.ndim(f)==3:
            f = transform.resize(f, [n, n, f.shape[2]], 1)
    return f

def perform_wavortho_transf(f, Jmin, dir, h):
    """
        perform_wavortho_transf - compute orthogonal wavelet transform

        fw = perform_wavortho_transf(f,Jmin,dir,options);

        You can give the filter in options.h.

        Works in 2D only.

        Copyright (c) 2014 Gabriel Peyre
    """

    n = f.shape[1]
    Jmax = int(np.log2(n)) - 1
    # compute g filter
    u = np.power(-np.ones(len(h) - 1), range(1, len(h)))
    # alternate +1/-1
    g = np.concatenate(([0], h[-1:0:-1] * u))

    if dir == 1:
        ### FORWARD ###
        fW = f.copy()
        for j in np.arange(Jmax, Jmin - 1, -1):
            A = fW[:2 ** (j + 1):, :2 ** (j + 1):]
            for d in np.arange(1, 3):
                Coarse = subsampling(cconv(A, h, d), d)
                Detail = subsampling(cconv(A, g, d), d)
                A = np.concatenate((Coarse, Detail), axis=d - 1)
            fW[:2 ** (j + 1):, :2 ** (j + 1):] = A
        return fW
    else:
        ### BACKWARD ###
        fW = f.copy()
        f1 = fW.copy()
        for j in np.arange(Jmin, Jmax + 1):
            A = f1[:2 ** (j + 1):, :2 ** (j + 1):]
            for d in np.arange(1, 3):
                if d == 1:
                    Coarse = A[:2**j:, :]
                    Detail = A[2**j: 2**(j + 1):, :]
                else:
                    Coarse = A[:, :2 ** j:]
                    Detail = A[:, 2 ** j:2 ** (j + 1):]
                Coarse = cconv(upsampling(Coarse, d), nt.reverse(h), d)
                Detail = cconv(upsampling(Detail, d), nt.reverse(g), d)
                A = Coarse + Detail
            f1[:2 ** (j + 1):, :2 ** (j + 1):] = A
        return f1


def plot_wavelet(fW, Jmin=0):
    """
        plot_wavelet - plot wavelets coefficients.

        U = plot_wavelet(fW, Jmin):

        Copyright (c) 2014 Gabriel Peyre
    """
    def rescaleWav(A):
        v = abs(A).max()
        B = A.copy()
        if v > 0:
            B = .5 + .5 * A / v
        return B
    ##
    n = fW.shape[1]
    Jmax = int(np.log2(n)) - 1
    U = fW.copy()
    for j in np.arange(Jmax, Jmin - 1, -1):
        U[:2 ** j:,    2 ** j:2 **
            (j + 1):] = rescaleWav(U[:2 ** j:, 2 ** j:2 ** (j + 1):])
        U[2 ** j:2 ** (j + 1):, :2 **
          j:] = rescaleWav(U[2 ** j:2 ** (j + 1):, :2 ** j:])
        U[2 ** j:2 ** (j + 1):, 2 ** j:2 ** (j + 1):] = (
            rescaleWav(U[2 ** j:2 ** (j + 1):, 2 ** j:2 ** (j + 1):]))
    # coarse scale
    U[:2 ** j:, :2 ** j:] = nt.rescale(U[:2 ** j:, :2 ** j:])
    # plot underlying image
    imageplot(U)
    # display crosses
    for j in np.arange(Jmax, Jmin - 1, -1):
        plt.plot([0, 2 ** (j + 1)], [2 ** j, 2 ** j], 'r')
        plt.plot([2 ** j, 2 ** j], [0, 2 ** (j + 1)], 'r')
    # display box
    plt.plot([0, n], [0, 0], 'r')
    plt.plot([0, n], [n, n], 'r')
    plt.plot([0, 0], [0, n], 'r')
    plt.plot([n, n], [0, n], 'r')
    return U

def psnr(x, y, vmax=-1):
    """
     psnr - compute the Peack Signal to Noise Ratio

       p = psnr(x,y,vmax);

       defined by :
           p = 10*log10( vmax^2 / |x-y|^2 )
       |x-y|^2 = mean( (x(:)-y(:)).^2 )
       if vmax is ommited, then
           vmax = max(max(x(:)),max(y(:)))

       Copyright (c) 2014 Gabriel Peyre
    """

    if vmax < 0:
        m1 = abs(x).max()
        m2 = abs(y).max()
        vmax = max(m1, m2)
    d = np.mean((x - y) ** 2)
    return 10 * np.log10(vmax ** 2 / d)

def snr(x, y):
    """
    snr - signal to noise ratio

       v = snr(x,y);

     v = 20*log10( norm(x(:)) / norm(x(:)-y(:)) )

       x is the original clean signal (reference).
       y is the denoised signal.

    Copyright (c) 2014 Gabriel Peyre
    """

    return 20 * np.log10(pylab.norm(x) / pylab.norm(x - y))

def subsampling(x, d):
    # subsampling along dimension d by factor p=2
    p = 2
    if d == 1:
        y = x[::p, :]
    elif d == 2:
        y = x[:, ::p]
    else:
        raise Exception('Not implemented')
    return y

def upsampling(x, d):
    """
        up-sampling along dimension d by factor p=2
    """
    p = 2
    s = x.shape
    if d == 1:
        y = np.zeros((p * s[0], s[1]))
        y[::p, :] = x
    elif d == 2:
        y = np.zeros((s[0], p * s[1]))
        y[:, ::p] = x
    else:
        raise Exception('Not implemented')
    return y


def plot_dictionary(D, title='Dictionary'):
    ''' Plot a dictionary of shape (width*width, n_atoms) '''
    # Check that D.shape == (width*width, n_atoms)
    assert len(D.shape) == 2
    assert int(np.sqrt(D.shape[0]))**2 == D.shape[0]
    (signal_size, n_atoms) = D.shape
    # Rescale values in each atom to have a max absolute value of 1
    # This gives brighter plots
    D = D / np.max(abs(D), axis=0)

    # Reshape dictionary to patches
    width = int(np.sqrt(D.shape[0]))
    D = D.reshape((width, width, n_atoms))
    n = int(np.ceil(np.sqrt(n_atoms)))  # Size of the plot in number of atoms

    # Pad the atoms
    pad_size = 1
    missing_atoms = n ** 2 - n_atoms

    padding = (((pad_size, pad_size), (pad_size, pad_size),
                (0, missing_atoms)))
    D = np.pad(D, padding, mode='constant', constant_values=1)
    padded_width = width + 2*pad_size
    D = D.reshape(padded_width, padded_width, n, n)
    D = D.transpose(2, 0, 3, 1)  # Needed for the reshape
    big_image_size = n*padded_width
    D = D.reshape(big_image_size, big_image_size)
    imageplot(D)
    plt.title(title)
    plt.show()
