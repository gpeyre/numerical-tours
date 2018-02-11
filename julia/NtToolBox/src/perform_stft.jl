function perform_stft(x,w,q,n)
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
    if ndims(x) == 1
        dir = 1
    else
        dir = -1
    end

    # perform sampling
    X = Array(1:q:n+2)

    p = length(X)
    eta = 1

    if w%2 == 1
        w = ceil((w-1)/2)*2+1
        w1 = Base.div(w-1,2)
        dX = Array(-w1:w1)
    else
        dX = Array((-Base.div(w,2)+1):Base.div(w,2))
    end

    X1 = repeat(X,outer=(1,w)) + repeat(dX, outer=(1,p))'
    #periodic boundary conditions
    X1 = ((X1-1).%n)+1;

    I = X1'
    I[I.<=0]+=n


    # build a sin weight function
    W = .5 *(1 - cos(2*pi*Array(0:(w-1))./(w-1)))

    #renormalize the windows
    weight = zeros(n)

    for i in 1:p
        weight[I[:,i]] = weight[I[:,i]] + W.^2
    end

    weight = sqrt(weight)
    Weight = repeat(W, outer=(1,p))

    for i in 1:p
        Weight[:,i] = Weight[:,i]./weight[I[:,i]]
    end

    #compute the transform
    if dir == 1
        y = zeros(eta*w,p)
        if w%2 == 1
            m = Base.div(eta*w+1,2)
            w1 = Base.div((w-1),2)
            sel = Array(m-w1:m+w1) - 1
        else
            m = Base.div((eta*w),2)+1
            w1 = Base.div(w,2)
            sel = Array(m-w1:m+w1-1) - 1
        sel[sel.<=0] += eta*w
        y[sel,:] = x[I].*Weight

        #perform the transform
        y = my_transform(y,+1)
        end
    else
        x = my_transform(x,-1)
        x = real.(x.*Weight)
        y = zeros(n)
        for i in 1:p
            y[I[:,i]] = y[I[:,i]] + x[:,i]
        end
    end
    return y
end
########################################################################
function my_transform(x,dir)

    # my_transform - perform either FFT with energy conservation.
    # Works on array of size (w,w,a,b) on the 2 first dimensions.
    w = size(x)[1]
    if dir == 1
        y = fft(x,1)/sqrt(w)
    else
        y = ifft(x*sqrt(w),1)
    end
    return y
end
