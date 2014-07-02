import numpy as np
import pylab
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
# try to not make use of transform.resize
from skimage import transform

def temp():
	return 2;

def subsampling(x,d):
    # subsampling along dimension d by factor p=2
    p = 2;
    if d==1: y = x[::p,:];
    elif d==2: y = x[:, ::p];
    else: raise Exception('Not implemented');
    return y;

def upsampling(x,d):
    # up-sampling along dimension d by factor p=2
    p = 2;
    s = x.shape;
    if d==1: 
        y = np.zeros((p*s[0],s[1]));
        y[::p,:] = x;
    elif d==2: 
        y = np.zeros((s[0],p*s[1]));
        y[:,::p] = x;
    else: raise Exception('Not implemented');
    return y;

def reverse(x): 
	return x[::-1];

def circshift(x,k):   
    # return concatenate((x[k::], x[0:k:]));
    return np.roll(x,-k, axis=0);

def cconv(x,h,d):
    # circular convolution along dimension d. 
    # h should be small and with odd size  
    if d==2:
        # apply to transposed matrix
        return np.transpose(cconv(np.transpose(x),h,1));
    y = np.zeros(x.shape);
    p = len(h);
    pc = (p-1)/2;
    for i in range(0,p):
        y = y + h[i]*circshift(x,i-pc);  
    return y;

def rescale(f):
	v = f.max()-f.min();
	g = (f-f.min()).copy();
	if v>0:
		g = g/v;
	return g;
	
def imageplot(f, str='', sbpt=[]):
	# use nearest neighbor interpolation
	if sbpt!=[]:
		plt.subplot(sbpt[0],sbpt[1],sbpt[2]);
	imgplot = plt.imshow(f, interpolation='nearest'); 
	imgplot.set_cmap('gray'); 
	pylab.axis('off');
	if str!='':
		plt.title(str);
	
def load_image(name, n=-1):
	f = plt.imread(name);
	# turn into normalized grayscale image
	f = rescale( np.sum(f, axis=2) );
	# change the size of the image
	if n>0:
		f = transform.resize(f, [n,n], 1); 
	return f;
	
def plot_wavelet(fW, Jmin=0):
	"""
		plot_wavelet - plot wavelets coefficients.
		
		U = plot_wavelet(fW, Jmin):
		
		Copyright (c) 2014 Gabriel Peyre
	"""
	def rescaleWav(A):
		v = abs(A).max();
		B = A.copy();
		if v>0:
			B = .5+.5*A/v;
		return B;
	n = fW.shape[1];
	Jmax = np.log2(n)-1;
	U = fW.copy();
	for j in np.arange(Jmax,Jmin-1,-1):
		U[    :2**j:,    2**j:2**(j+1):] = rescaleWav(U[:2**j:,2**j:2**(j+1):]);
		U[2**j:2**(j+1):,    :2**j:] = rescaleWav(U[2**j:2**(j+1):,:2**j:]);
		U[2**j:2**(j+1):,2**j:2**(j+1):] = rescaleWav(U[2**j:2**(j+1):,2**j:2**(j+1):]);
	# coarse scale 
	U[:2**j:,:2**j:] = rescale(U[:2**j:,:2**j:]);
	# plot underlying image
	imageplot(U);
	# display crosses
	for j in np.arange(Jmax,Jmin-1,-1):
		plt.plot([0, 2**(j+1)], [2**j, 2**j], 'r');
		plt.plot([2**j, 2**j], [0, 2**(j+1)], 'r');
	# display box
	plt.plot([0, n], [0, 0], 'r');
	plt.plot([0, n], [n, n], 'r');
	plt.plot([0, 0], [0, n], 'r');
	plt.plot([n, n], [0, n], 'r');
	return U;
	
def perform_wavortho_transf(f,Jmin,dir,h):
	"""" 
		perform_wavortho_transf - compute orthogonal wavelet transform
	
	   	fw = perform_wavortho_transf(f,Jmin,dir,options);
	
	   	You can give the filter in options.h.
	
	   	Works in 2D only.
	
	   	Copyright (c) 2014 Gabriel Peyre
	"""
	
	n = f.shape[1];
	Jmax = np.log2(n)-1;	
	# compute g filter
	u = np.power(-np.ones(len(h)-1),range(1,len(h))); # alternate +1/-1
	g = np.concatenate(([0], h[-1:0:-1] * u));
	
	if dir==1:
	    ### FORWARD ###
		fW = f.copy();
		for j in np.arange(Jmax,Jmin-1,-1):
		    A = fW[:2**(j+1):,:2**(j+1):];
		    for d in np.arange(1,3):
		        Coarse = subsampling(cconv(A,h,d),d);
		        Detail = subsampling(cconv(A,g,d),d);
		        A = np.concatenate( (Coarse, Detail), axis=d-1 );
		    fW[:2**(j+1):,:2**(j+1):] = A;
		return fW;
	else:
		### BACKWARD ###
		fW = f.copy();
		f1 = fW.copy();
		for j in np.arange(Jmin,Jmax+1): 
		    A = f1[:2**(j+1):,:2**(j+1):];
		    for d in np.arange(1,3):
		        if d==1:
		            Coarse = A[:2**j:,:];
		            Detail = A[2**j:2**(j+1):,:];
		        else:
		            Coarse = A[:,:2**j:];
		            Detail = A[:,2**j:2**(j+1):];               
		        Coarse = cconv(upsampling(Coarse,d),reverse(h),d);
		        Detail = cconv(upsampling(Detail,d),reverse(g),d);
		        A = Coarse + Detail;
		    f1[:2**(j+1):,:2**(j+1):] = A;
		return f1;

def snr(x,y):
		"""
			snr - signal to noise ratio
		
		   v = snr(x,y);
		
		 v = 20*log10( norm(x(:)) / norm(x(:)-y(:)) )
		
		   x is the original clean signal (reference).
		   y is the denoised signal.
		 
		   Copyright (c) 2014 Gabriel Peyre
		"""

		return 20*np.log10(pylab.norm(x)/pylab.norm(x-y));
		
def psnr(x,y, vmax=-1):
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
		
		if vmax<0:
		    m1 = abs(x).max();
		    m2 = abs(y).max();
		    vmax = max(m1,m2);
		d = np.mean( (x-y)**2 );
		return 10*np.log10( vmax**2/d );

def clamp(x,a=[],b=[]):
		"""
		 clamp - clamp a value
		
		   y = clamp(x,a,b);
		
		 Default is [a,b]=[0,1].
		
		   Copyright (c) 2004 Gabriel Peyre
		"""
		
		if a==[]:
		    a = 0.0;
		if b==[]:
		    b = 1.0;
		return np.minimum(np.maximum(x,a),b);