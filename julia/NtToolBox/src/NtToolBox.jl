
module NtToolBox

# package code goes here

export load_image, imageplot, rescale, clamP, snr, perform_wavelet_transf, plot_wavelet, compute_wavelet_filter, perform_wavortho_transf, grad, div, plot_levelset, gaussian_blur, bilinear_interpolate, Grad, Div, perform_redistancing, perform_haar_transf, perform_fast_marching, load_sound, perform_stft, plot_spectrogram, compute_max, perform_blurring, plot_vf, meshgrid, read_mesh, compute_boundary, plot_mesh, compute_normal, perform_linprog, plot_hufftree, perform_conjugate_gradient, perform_dijkstra_fm

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

function compute_max(X,d)
    """ compute_max - compute maximum along dimension d

       Y,I = compute_max(X,d);

       Copyright (c) 2008 Gabriel Peyre
    """
    if ndims(X)>=2
        Y = maximum(X,d)
        I = mapslices(indmax,X,d)[:]
    elseif ndims(X)==1
        Y = maximum(X)
        I = indmaw(X)
    end
    return Y,I
end

function compute_max(X)
    if ndims(X)>=2
        Y = maximum(X,1)
        I = mapslices(indmax,X,1)[:]
    elseif ndims(X)==1
        Y = maximum(X)
        I = indmaw(X)
    end
    return Y,I
end


# function circshift(x, p)
#     """
#         Circular shift of an array.
#     """
#     y = copy(x)
#     y = cat(1, y[p[1]+1:end, :], y[1:p[1]+1, :])
#     if (size(x)[2] > 0) & (length(p) > 1)
#         y = cat(2, y[:, p[1]+1:end], y[:, 1:p[1]+1])
#     end
#     return y
# end


include("ndgrid.jl")
include("signal.jl")
include("perform_wavelet_transf.jl")
include("compute_wavelet_filter.jl")
include("Grad.jl")
include("Div.jl")
include("perform_blurring.jl")
include("read_bin.jl")
include("isosurface.jl")
include("perform_thresholding.jl")
include("load_sound.jl")
include("perform_stft.jl")
include("plot_spectrogram.jl")
include("plot_vf.jl")
include("read_mesh.jl")
include("compute_boundary.jl")
include("plot_mesh.jl")
include("compute_normal.jl")
include("perform_linprog.jl")
include("plot_hufftree.jl")
include("perform_conjugate_gradient.jl")
include("graph.jl")
include("perform_redistancing.jl")
include("perform_haar_transf.jl")
include("perform_fast_marching.jl")



end # module
