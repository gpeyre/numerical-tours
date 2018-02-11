w = randn(n0,n0)/n0
w = w-mean(vec(w))+1/n0^2
f = real(ifft(fft(f1).*repeat(fft(w), outer=(1, 1, size(f0,3)))))
clf()
u = ones(n0,n0,3)
range = Base.div(n0,2)-Base.div(n,2)+1:Base.div(n0,2) + Base.div(n,2)
u[range, range, :] = f0
u[1,1,:]=0; u[2,1,:]=1
imageplot(clamp(u,0,1), "Input", [1,2,1]);
u = f; u[1,1,:]=0; u[2,1,:]=1
imageplot(clamp(u,0,1), "Synthesized", [1,2,2]);
