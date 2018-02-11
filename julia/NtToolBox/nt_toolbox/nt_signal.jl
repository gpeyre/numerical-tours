module nt_signal

using PyPlot

export imageplot, resize_image, load_image, cconvol, upsampling, subsampling, plot_wavelet, snr, perform_wavortho_transf

#######################################
function imageplot(g, str="", a=-1, b=-1, c=-1)
	if size(g,3)==1
		# special for b&w images.
		g = g[:,:,1];
	end
	if c>0
		subplot(a,b,c);
	end
	imshow(g, interpolation="nearest")
	set_cmap("gray")
	axis("off")
	if str!=""
		title(str);
	end
end

#######################################
function load_image(name,N)
	g = imread(name);
	g = resize_image(g, N);
	return g;
end

#######################################
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
function plot_wavelet(fW, Jmin=0)
    #    plot_wavelet - plot wavelets coefficients.
    #    U = plot_wavelet(fW, Jmin):
    #    Copyright (c) 2014 Gabriel Peyre

    function rescaleWav(A)
        v = maximum( abs(A[:]) );
        B = copy(A)
        if v > 0
            B = .5 + .5 * A / v;
        end
        return B;
    end
    function rescale(x,a=0,b=1)
		m = minimum(x[:]);
		M = maximum(x[:]);
		if M-m<1e-10
		    y = x;
		else
		    y = (b-a) * (x-m)/(M-m) + a;
		end
		return y;
	end ## end of rescale
    ## 
    n = size(fW,1);
    Jmax = log2(n) - 1;
    U = copy(fW);
    for j = Jmax:-1:Jmin
    	L = 1:2^j; H = (2^j)+1:2^(j+1);
        U[L,H] = rescaleWav( U[L,H] );
        U[H,L] = rescaleWav( U[H,L] );
        U[H,H] = rescaleWav( U[H,H] );
    end
    # coarse scale
    L = 1:2^Jmin;
    U[L,L] = rescale(U[L,L]);
    # plot underlying image
    imageplot(U);
    # display crosses
    for j=Jmax:-1:Jmin
        plot([0, 2^(j+1)], [2^j, 2^j], "r");
        plot([2^j, 2^j], [0, 2^(j + 1)], "r");
    end
    # display box
    plot([0, n], [0, 0], "r");
    plot([0, n], [n, n], "r");
    plot([0, 0], [0, n], "r");
    plot([n, n], [0, n], "r");
    return U;

end # plot_wavelet



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
	Jmax = log2(n)-1; 

	g = [0, h[length(h):-1:2]] .* (-1).^(1:length(h));

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

#####################################################
function snr(x,y)

	# snr - signal to noise ratio
	#
	#   v = snr(x,y);
	#
	# v = 20*log10( norm(x[:]) / norm(x[:]-y[:]) )
	#
	#   x is the original clean signal (reference).
	#   y is the denoised signal.
	# 
	#   Copyright (c) 2008 Gabriel Peyre

	return 20*log10(norm(x[:])/norm(x[:]-y[:]));

end # snr

end ## module nt_signal
