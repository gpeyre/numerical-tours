
module NtToolBox

# package code goes here

export load_image, imageplot, rescale, clamP, snr, perform_wavelet_transf, plot_wavelet, compute_wavelet_filter, perform_wavortho_transf, grad, div, plot_levelset, gaussian_blur, bilinear_interpolate, Grad, Div

## m must be of type Array{Float32, 3} or Array{Float64, 3}, it makes products between the vectors of the 3rd dimension and v





function Mdot(m, v::Array{Float64, 1})
  w = ones(size(m, 1), size(m, 2))
  for i in 1:size(m, 1)
    for j in 1:size(m, 2)
      w[i, j] = (m[i, j, :]'*v)[1]
    end
  end
  return w
end

##

## Rescale linearly the dynamic of a vector to fit within a range [a,b]

function rescale(f, a = 0, b = 1)
  v = maximum(f) - minimum(f)
  g = copy(f - minimum(f))
  if v > 0
    g = g / v
  end
  return a + g*(b - a)
end

##

## Function clamp

function clamP(x, a = [], b = [])

  """
   clamp - clamp a value

     y = clamp(x,a,b);

   Default is [a,b]=[0,1].

     Copyright (c) 2004 Gabriel Peyre
  """

  if a == []
    a = 0.0
  end
  if b == []
    b = 1.0
  end
  return min(max(x, a), b)
end

#function plot_levelset(Z, level = 0, f = [])
#    """
#        f is supposed to be of the same shape as Z
#    """
#     if length(f) == 0
#         f = copy(Z)
#     end
#
#     (n, p) = size(Z)
#     (X, Y) = meshgrid(collect(0 : n - 1), collect(0 : p - 1))
#     contour(X, Y, Z, [level], linewidths = 2, colors = "red")
#     imageplot(f)
# end

function gaussian_blur(f, sigma)

      """ gaussian_blur - gaussian blurs an image
      %
      %   M = perform_blurring(M, sigma, options);
      %
      %   M is the original data
      %   sigma is the std of the Gaussian blur (in pixels)
      %
      %   Copyright (c) 2007 Gabriel Peyre
      """
      if sigma <= 0
          return
      end
      n = maximum(size(f))
      t = [collect(0:n/2); collect(-n/2:-2)]
      (Y, X) = meshgrid(t, t)
      h = exp( -(X.^2 + Y.^2)./(2.0*float(sigma)^2) )
      h = h./sum(h)
      return real(plan_ifft((plan_fft(f)*f).*(plan_fft(h)*h))*((plan_fft(f)*f).*(plan_fft(h)*h)))
end

function plot_levelset(Z, level = 0, f = [])
    """
        f is supposed to be of the same shape as Z
    """
    if length(f) == 0
        f = copy(Z)
    end

    (n, p) = size(Z)
    (X, Y) = meshgrid(collect(0 : n - 1), collect(0 : p - 1))
    contour(X, Y, Z, [level], linewidths = 2, colors="red")
    imageplot(f)
end

function Grad(M, bound = "sym", order = 1)
    """
        grad - gradient, forward differences

          [gx,gy] = grad(M, options);
        or
          g = grad(M, options);

          options.bound = 'per' or 'sym'
          options.order = 1 (backward differences)
                        = 2 (centered differences)

          Works also for 3D array.
          Assme that the function is evenly sampled with sampling step 1.

          See also: div.

          Copyright (c) Gabriel Peyre
    """


    # retrieve number of dimensions
    nbdims = ndims(M)


    if bound == "sym"
        nx = size(M)[1]
        if order == 1
            fx = M[hcat((collect(2 : nx),[nx])), :] - M
        else
            fx = (M[hcat((collect(2 : nx), [nx])), :] - M[hcat(([1], collect(1 : nx - 1))), :])./2.
            # boundary
            fx[1, :] = M[2, :] - M[1, :]
            fx[nx, :] = M[nx, :] - M[nx - 1, :]
        end

        if nbdims >= 2
            ny = size(M)[2]
            if order == 1
                fy = M[:, hcat((collect(2 : ny), [ny]))] - M
            else
                fy = (M[:, hcat((collect(2 : ny), [ny]))] - M[:, hcat(([1], collect(1 : ny - 1)))])./2.
                # boundary
                fy[:, 1] = M[:, 2] - M[:, 1]
                fy[:, ny] = M[:, ny]-M[:, ny - 1]
            end
        end

        if nbdims >= 3
            nz = size(M)[3]
            if order == 1
                fz = M[:, :, hcat((collect(2 : nz), [nz]))] - M
            else
                fz = (M[:, :, hcat((collect(2 : nz), [nz]))] - M[:, :, hcat(([1], collect(1 : nz - 1)))])./2.
                # boundary
                fz[:, :, 1] = M[:, :, 2] - M[:, :, 1]
                fz[:, :, ny] = M[:, :, nz] - M[:, :, nz - 1]
            end
        end
    else
        nx = size(M)[1]
        if order == 1
            fx = M[hcat((collect(2 : nx), [1])), :] - M
        else
            fx = (M[hcat((collect(2 : nx), [1])), :] - M[hcat(([nx], collect(1 : nx - 1))), :])./2.
        end

        if nbdims >= 2
            ny = size(M)[2]
            if order == 1
                fy = M[:, hcat((collect(2 : ny), [1]))] - M
            else
                fy = (M[:, hcat((collect(2 : ny), [1]))] - M[:, hcat(([ny], collect(1 : ny - 1)))])./2.
            end
        end

        if nbdims >= 3
            nz = size(M)[3]
            if order == 1
                fz = M[:, :, hcat((collect(2 : nz), [1]))] - M
            else
                fz = (M[:, :, hcat((collect(2 : nz), [1]))] - M[:, :, hcat(([nz], collect(1 : nz - 1)))])./2.
            end
        end
    end

    if nbdims==2
        fx = cat(3, fx[:, :], fy[:, :])
    elseif nbdims==3
        fx = cat(4, (fx[:, :, :], fy[:, :, :], fz[:, :, :]))
    end

    return fx
end


include("ndgrid.jl")
include("signal.jl")
include("perform_wavelet_transf.jl")
include("compute_wavelet_filter.jl")
include("Div.jl")
include("perform_blurring.jl")




end # module
