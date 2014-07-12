f = real(ifft2(fft2(f0).*repmat(fft2(w), [1 1 size(f0,3)])));
clf;
u = f0; u(1,1,:)=0; u(2,1,:)=1;
imageplot(clamp(u), 'Input', 1,2,1);
u = f; u(1,1,:)=0; u(2,1,:)=1;
imageplot(clamp(u), 'Synthesized', 1,2,2);
