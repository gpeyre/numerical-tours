using PyPlot

## Rescale linearly the dynamic of a vector to fit within a range [a,b]

function rescale(f, a = 0, b = 1)
  v = maximum(f) - minimum(f)
  g = copy(f - minimum(f))
  if v > 0
    g = g / v
  end
  return a + g*(b - a)
end    
