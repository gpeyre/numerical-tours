using PyPlot

## Rescale linearly the dynamic of a vector to fit within a range [a,b]

function rescale(f, a = 0, b = 1
  v = maximum(f) - minimum(f)
  g = copy(f - minimum(f))
  if v > 0
    g = g / v
  end
  return a + g*(b - a)
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
