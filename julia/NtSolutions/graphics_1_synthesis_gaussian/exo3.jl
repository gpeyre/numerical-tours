f = real(ifft(fft(f0).*repeat(fft(w), outer=(1, 1,size(f0,3)))));
u = f0; u[1,1,:]=0; u[2,1,:]=1
imageplot(clamp(u,0,1), "Input", [1,2,1])
u = f; u[1,1,:]=0; u[2,1,:]=1
imageplot(clamp(u,0,1), "Synthesized", [1,2,2]);
