using PyPlot
using Images
using NtToolBox.rescale
using NtToolBox.Mdot

##

## Load an image from a file, rescale its dynamic to [0,1], turn it into a grayscale image and resize it to size n x n.

function load_image(name, n = -1, flatten = 1, resc = 1, grayscale = 1)

  if ndims(PyPlot.imread(name)) >=3
    f = PyPlot.imread(name)[:, :, 1:3]
  end

  f = PyPlot.imread(name)

  if grayscale == 1
    if (flatten == 1) & (ndims(f) > 2)
      f = sum(f, 3)   # We must be careful, this returns an array of dimension 512*512*1 while
      f = f[:, :]     # in the python file np.sum(f, axis=2) returns an array of dimension 512*512

    end
  end
  if resc == 1
    f = rescale(f)  # rescale must be defined in the file general.jl of the package NtToolBox
                    # and we must find a way to import it to this file signal.jl
  end
  if n > 0
    if ndims(f) == 3
        f = Images.imresize(f,(n,n,3))
    elseif ndims(f) == 2
      f = Images.imresize(f, (n, n))
    end
  end
  return f
end

##

## Use nearest neighbor interpolation for the display.

function imageplot(f, str = "", sbpt = [])    # Function to check.

  if sbpt != []
    subplot(sbpt[1], sbpt[2], sbpt[3])
  end
  imshow(f, interpolation = "nearest")
  set_cmap("gray")
  axis("off")
  if str != ""
    title(str)
  end
end

function snr(x, y)
  """
  snr - signal to noise ratio

     v = snr(x,y);

   v = 20*log10( norm(x(:)) / norm(x(:)-y(:)) )

     x is the original clean signal (reference).
     y is the denoised signal.

  Copyright (c) 2014 Gabriel Peyre
  """
  return 20 * log10(norm(x[:]) / norm(x[:] - y[:]))
end

############################################################################

function resize_image(f, N)
	P = size(f,1);
	# add 1 pixel boundary
	g = f;
	g = cat(2, g, g[:,1,:]);
	g = cat(1, g, g[1,:,:]);
	# interpolate
	t = linspace(1,P,N);
	ti = floor(t); tj = ceil(t);
	fi = t-floor(t); fj = 1-fi;
	h = zeros(N,N,size(f,3));
	for s=1:size(f,3)
	    h[:,:,s] = g[ti,ti,s] .* (fj*fj') + g[tj,tj,s] .* (fi*fi') + g[ti,tj,s] .* (fj*fi') + g[tj,ti,s] .* (fi*fj');
	end
	return h;
end

###################################################
function trim_dim(x)
	# remove empty trailing dimensions
	if size(x,3)==1
		x = x[:,:,1];
		if size(x,2)==1
			x = x[:,1];
		end
	end
	return x;
end

###################################################
function subsampling(x,d=1,p=2)

	# downsampling - subsampling along dimension d
	#
	#   y = subsampling(x,d,p);
	#
	#   default is p==2, d==1
	#
	#   Copyright (c) 2009 Gabriel Peyre

	if d==1
	        y = x[1:p:end,:,:];
	elseif d==2
	        y = x[:,1:p:end,:];
	elseif d==3
	        y = x[:,:,1:p:end];
	else
	        error("Not implemented");
	end
	y = trim_dim(y);

	return y;

end # subsampling

###################################################
function upsampling(x,d=1,p=2)

	# upsampling - add p zeros between samples along dimension d
	#
	#   y = upsampling(x,d,p);
	#
	#   default is p==2, d==1
	#
	#   Copyright (c) 2009 Gabriel Peyre

	if d==1
	        y = zeros(p*size(x,1),size(x,2),size(x,3));
	        y[1:p:end,:,:] = x;
	elseif d==2
	        y = zeros(size(x,1),p*size(x,2),size(x,3));
	        y[:,1:p:end,:] = x;
	elseif d==3
	        y = zeros(size(x,1),size(x,2),p*size(x,3));
	        y[:,:,1:p:end] = x;
	else
	        error("Not implemented");
	end
	y = trim_dim(y);

	return y;
end # upsampling


###################################################
function cconvol(x,h,d=1)

	# cconvol - spatial domain circular convolution
	#
	#   y = cconv(x,h,d);
	#
	#   Circular convolution along dimension |d|.
	#
	#   If p=length(h) is odd then h((p+1)/2) corresponds to the zero position.
	#   If p is even, h(p/2) corresponds to the zero position.
	#
	#   Copyright (c) 2009 Gabriel Peyre

	if d==2
	    return ( cconvol( x',h, 1 )' );
	elseif d==3
		# permute do not exists in python
	    return permute( cconvol( permute(x,[3, 2, 1]),h), [3 2 1]);
	end

	p = length(h);
	if mod(p,2)==0
	    pc = p/2;
	else
	    pc = (p+1)/2;
	end

	y = zeros(size(x));
	for i=1:length(h)
	    y = y + h[i]*circshift(x,i-pc);
	end
	return y;

end # cconvol

#####################################################
# function plot_wavelet(fW, Jmin=0)
#     #    plot_wavelet - plot wavelets coefficients.
#     #    U = plot_wavelet(fW, Jmin):
#     #    Copyright (c) 2014 Gabriel Peyre
#
#     function rescaleWav(A)
#         v = maximum( abs(A[:]) );
#         B = copy(A)
#         if v > 0
#             B = .5 + .5 * A / v;
#         end
#         return B;
#     end
#     function rescale(x,a=0,b=1)
# 		m = minimum(x[:]);
# 		M = maximum(x[:]);
# 		if M-m<1e-10
# 		    y = x;
# 		else
# 		    y = (b-a) * (x-m)/(M-m) + a;
# 		end
# 		return y;
# 	  end ## end of rescale
#     ##
#     n = size(fW,1);
#     Jmax = Int(log2(n) - 1);
#     U = copy(fW);
#     for j = Jmax:-1:Jmin
#     	L = Array{Int64,1}(1:2^j); H = Array{Int64,1}((2^j)+1:2^(j+1));
#         U[L,H] = rescaleWav( U[L,H] );
#         U[H,L] = rescaleWav( U[H,L] );
#         U[H,H] = rescaleWav( U[H,H] );
#     end
#     # coarse scale
#     L = Array{Int64,1}(1:2^Jmin);
#     U[L,L] = rescale(U[L,L]);
#     # plot underlying image
#     imageplot(U);
#     # display crosses
#     for j=Jmax:-1:Jmin
#         plot([0, 2^(j+1)], [2^j, 2^j], "r");
#         plot([2^j, 2^j], [0, 2^(j + 1)], "r");
#     end
#     # display box
#     plot([0, n], [0, 0], "r");
#     plot([0, n], [n, n], "r");
#     plot([0, 0], [0, n], "r");
#     plot([n, n], [0, n], "r");
#     return U;
#
# end # plot_wavelet

#########
function plot_wavelet(fW, Jmin = 0)
    """
        plot_wavelet - plot wavelets coefficients.

        U = plot_wavelet(fW, Jmin):

        Copyright (c) 2014 Gabriel Peyre
    """
    function rescaleWav(A)
        v = maximum(abs(A[:]))
        B = copy(A)
        if v > 0
            B = .5 + .5 .* A ./ v
        end
        return B
    end
    ##
    n = size(fW)[2]
    Jmax = Int(log2(n) - 1)
    U = copy(fW)
    for j in Jmax:-1:Jmin
        U[1:2^j,    2^j+1:2 ^
            (j + 1)] = rescaleWav(U[1:2^j, 2^j+1:2^(j + 1)])
        U[2^j+1:2^(j+1), 1:2 ^
          j] = rescaleWav(U[2^j+1:2^(j+1), 1:2^j])
        U[2^j+1:2^(j+1), 2^j+1:2^(j+1)] = (
            rescaleWav(U[2^j+1:2^(j+1), 2^j+1:2^(j+1)]))
    end

    # coarse scale
    U[1:2^Jmin, 1:2^Jmin] = rescale(U[1:2^Jmin, 1:2^Jmin])

    # plot underlying image
    imageplot(U)
    # display crosses
    for j in Jmax:-1:Jmin
        plot([0, 2^(j+1)], [2^j, 2^j], "r")
        plot([2^j, 2^j], [0, 2^(j+1)], "r")
    end
    # display box
    plot([0, n], [0, 0], "r")
    plot([0, n], [n, n], "r")
    plot([0, 0], [0, n], "r")
    plot([n, n], [0, n], "r")
    return U
end

########



#########################################
function subselectdim(f,sel,d)
	g = [];
	if d==1
	    g = f[sel,:,:];
	elseif d==2
	    g = f[:,sel,:];
	elseif d==3
	    g = f[:,:,sel];
	end
	return trim_dim(g);
end


#########################################
function subassign(f,sel,g)
	d = ndims(f);
	if d==1
	    f[sel] = g;
	elseif d==2
	    f[sel,sel] = g;
	elseif d==3
	    f[sel,sel,sel] = g;
	end
	return f;
end # subassign


#########################################
function subselect(f,sel)
	d = ndims(f)
    if d==1
        return f[sel];
    elseif d==2
        return f[sel,sel];
    elseif d==3
        return f[sel,sel,sel];
	end
	return [];
end


##############################################################
function perform_wavortho_transf(f0,Jmin,dir,h)

	# perform_wavortho_transf - compute orthogonal wavelet transform
	#
	#   fw = perform_wavortho_transf[f,Jmin,dir,options);
	#
	#   You can give the filter in options.h.
	#
	#   Works in arbitrary dimension.
	#
	#   Copyright (c) 2009 Gabriel Peyre

	n = size(f0,1);
	Jmax = Int(log2(n)-1);

	g = [0; h[length(h):-1:2]] .* (-1).^(1:length(h)); # Je pense qu'il faut changer [0, h[length(h):-1:2]] en [0; h[length(h):-1:2]]

	f = copy(f0)
	if dir==1
	    ### FORWARD ###
	    for j=Jmax:-1:Jmin
	        sel = 1:2^(j+1);
	        a = subselect(f,sel);
	        for d=1:ndims(f)
	            a = cat(d, subsampling(cconvol(a,h,d),d), subsampling(cconvol(a,g,d),d) );
	        end
	        f = subassign(f,sel,a);
	    end
	else
	    ### FORWARD ###
	    for j=Jmin:Jmax
	        sel = 1:2^(j+1);
	        a = subselect(f,sel);
	        for d=1:ndims(f)
	            w = subselectdim(a,2^j+1:2^(j+1),d);
	            a = subselectdim(a,1:2^j,d);
	            a = cconvol(upsampling(a,d),reverse(h),d) + cconvol(upsampling(w,d),reverse(g),d);
	        end
	        f = subassign(f,sel,a);
	    end
	end
	return f;
end # perform_wavortho_transf

#########################################################

function div(g)
  """
      Compute a finite difference approximation of the gradient of a 2D vector field, assuming periodic BC.
  """

  S = size(g)
  s0 = [[S[1]]; collect(1:S[1] - 1)]
  s1 = [[S[2]]; collect(1:S[2] - 1)]
  f = (g[:, :, 1] - g[s0, :, 1]) + (g[:, : , 2] - g[:, s1, 2])
  return f
end

function grad(f)
    """
        Compute a finite difference approximation of the gradient of a 2D image, assuming periodic BC.
    """
    S = size(f)
    s0 = [collect(2:S[1]); [1]]
    s1 = [collect(2:S[2]); [1]]
    g = cat(3, f[s0, :] - f, f[:, s1] - f)
    return g
end

function bilinear_interpolate(im, x, y)

    # x0 = Array{Int64,1}(floor(x))
    x0 = [Int(floor(x))]
    x1 = x0 + 1
    # y0 = Array{Int64,1}(floor(y))
    y0 = [Int(floor(y))]
    y1 = y0 + 1

    x0 = clamp(x0, 1, size(im)[2])
    x1 = clamp(x1, 1, size(im)[2])
    y0 = clamp(y0, 1, size(im)[1])
    y1 = clamp(y1, 1, size(im)[1])

    Ia = im[ y0, x0 ]
    Ib = im[ y1, x0 ]
    Ic = im[ y0, x1 ]
    Id = im[ y1, x1 ]

    wa = (x1 - x) .* (y1 - y)
    wb = (x1 - x) .* (y - y0)
    wc = (x - x0) .* (y1 - y)
    wd = (x - x0) .* (y - y0)

    return wa.*Ia + wb.*Ib + wc.*Ic + wd.*Id
end
