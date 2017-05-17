clf()
for i in 1:4
    w = randn(n,n)/n
    w = w-mean(vec(w))+1/n^2
    f = real(ifft(fft(f0).*repeat(fft(w), outer=(1, 1, size(f0,3)))))
    u = f; u[1,1,:]=0; u[2,1,:]=1
    subplot(2,2,i)
    imageplot(clamp(u,0,1))
end
