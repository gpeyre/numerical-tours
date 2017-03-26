cR = sort(reshape(abs(fF), size(fF)[1]*size(fF)[2]))[end : -1 : 1]
h = plot(log10(cR), linewidth = 2)
xlim(0, n^2)

show()
