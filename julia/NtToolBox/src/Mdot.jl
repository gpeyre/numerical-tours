function Mdot(m, v::Array{Float64, 1})
  w = ones(size(m, 1), size(m, 2))
  for i in 1:size(m, 1)
    for j in 1:size(m, 2)
      w[i, j] = (m[i, j, :]'*v)[1]
    end
  end
  return w
end
