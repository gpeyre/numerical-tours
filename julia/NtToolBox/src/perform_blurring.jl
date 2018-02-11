
function perform_blurring(M, sigma, bound="sym")
    """
        perform_blurring - gaussian blurs an image

        M = perform_blurring(M, sigma, options);

        M is the original data
        sigma is the width of blurs (in pixels)

        Copyright (c) 2007 Gabriel Peyre
    """

    if sigma == zeros(size(sigma))
        return M
    end

    if ndims(M) > 2
        for i in 1 : size(M)[3]
             M[:, :, i] = perform_blurring(M[:, :, i], sigma, bound)
        end
    end

    n = maximum(size(M))

    eta = 4
    p = round((sigma.*eta)./2.).*2 + 1
    p = min(p, (round(n/2.)*2 - 1).*ones(length(p)))

    A = [1., 1.]
    if ndims(M) == 1
        A = 1 #1D
    end

    h = compute_gaussian_filter(p.*A, sigma./(4.*n), n*A)
    M = perform_convolution(M, h, bound)

    return M
end


function compute_gaussian_filter(n, s, N)
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
    nd = 1
    if (length(n) > 1) & (n[2] > 1)
        nd = 2
    end

    if (nd == 2) & (length(s) == 1)
        s = vcat(s, s)
    end

    if (nd == 2) & (length(N) == 1)
        N = vcat(N, N)
    end

    if nd == 1
        f = build_gaussian_filter_1d(n, s, N)
    else
        f = build_gaussian_filter_2d(n, s, N)
    end
    return f
end

function build_gaussian_filter_2d(n, s, N = [])
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

    # n = np.asarray(n)
    # s = np.asarray(s)
    # N = np.asarray(N)

    if length(N) == 0
        N = n
    end

    if (length(N) == 1) | (N[1] == 1)
        N = vcat(N, N)
    end

    if (length(s) == 1) | (s[1] == 1)
        s = vcat(s,s)
    end

    if length(collect(x for x in s if x <= 0)) > 0
        f = zeros(n)
        f[Int(round((n - 1)/2))] = 1
        return f
    end

    x = (collect(0 : n[1]) .- (n[1] - 1)/2.)/(N[1] - 1)
    y = (collect(0 : n[2]) .- (n[2] - 1)/2.)/(N[2] - 1)
    (Y, X) = meshgrid(y, x) #meshgrid is in the file ndgrid.jl
    f = exp(-(X.^2 ./ (2*s[1]^2)) - (Y.^2 ./ (2*s[2]^2)))
    f = f./sum(f)
    return f
end


function build_gaussian_filter_1d(n, s, N=[])
    """
        build_gaussian_filter_1d - compute a Gaussian filter.

        f = build_gaussian_filter_1d(n,s,N);

        Copyright (c) 2004 Gabriel Peyre
    """
    if length(N) == 0
        N = n
    end

    n = n[1]
    s = s[1]
    N = N[1]

    if s <= 0
        f = zeros(n)
        f[round((n - 1)/2)] = 1
        return f
    end

    x = (collect(0 : n) .- (n - 1)/2.)/(N - 1)
    f = exp(-x.^2 ./ (2*s^2))
    f = f./sum(f)
    return f
end


function perform_convolution(x, h, bound = "sym")
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

    if !(bound in ["sym", "per"])
        error("bound should be sym or per")
    end

    if ndims(x) == 3
      if size(x)[3] < 4
          #for color images
        y = x
        for i in 1 : size(x)[3]
          y[:, :, i] = perform_convolution(x[:, :, i], h, bound)
        end
      end
    return y
    end

    if ndims(x) == 3
      if size(x)[3] >= 4
        error("Not yet implemented for 3D array, use smooth3 instead.")
      end
    end

    n = size(x)
    p = size(h)

    nd = ndims(x)

    if nd == 1
        n = length(x)
        p = length(h)
    end

    if bound == "sym"

                #################################
        # symmetric boundary conditions #
        d1 = Array{Int64}(p)./2  # padding before
        d2 = p - d1 - 1    			    # padding after

        if nd == 1
        ################################# 1D #################################
            nx = length(x)
            xx = hcat(x[d1: -1: end], x, x[nx - 1 : -1 : nx - d2 - 1]) #A vérifier
            y = conv2(xx,h)
            y = y[p + 1 : nx - p - 1] #Il faut revérifier les indices

        elseif nd == 2
        ################################# 2D #################################
            #double symmetry
            nx = size(x)
            ny = size(x)
            xx = x
            xx = hcat(xx[d1[1] + 1 : -1 : end, :], xx, xx[nx : -1 : nx - d2[2] - 1, :])
            xx = vcat(xx[:, d1[2] : -1 : end], xx, xx[:, ny : -1 : ny - d2[2] - 1])
            y = conv2(xx, h) #mode="same"
            y = y[(2*d1[1] + 1) : (2*d1[1] + n[1] + 1), (2*d1[2] + 1) : (2*d1[2] + n[2] + 1)]
        end

    else

        ################################
        # periodic boundary conditions #

        if p > n
            error("h filter should be shorter than x.")
        end
        n = [n[1], n[2]]
        p = [p[1], p[2]]
        d = Array{Int64, 1}(floor((p - 1)/2.))
        if nd == 1
            h = hcat(h[d : end], hcat(zeros(n - p), h[1 : d - 1]))
            y = real(ifft((fft(x)*x).*(fft(h)*h)*((fft(x)*x).*(fft(h)*h))))
        else
            h = vcat(h[d[1] + 1 : end, :], vcat(zeros(n[1] - p[1], p[2]), h[1 : d[1], :]))
            h = hcat(h[:, d[2] + 1 : end], hcat(zeros(n[1], n[2] - p[2]), h[:, 1 : d[2] ]))
            y = real(plan_ifft((plan_fft(x)*x).*(plan_fft(h)*h))*((plan_fft(x)*x).*(plan_fft(h)*h)))
        end
    end
    return y
end
