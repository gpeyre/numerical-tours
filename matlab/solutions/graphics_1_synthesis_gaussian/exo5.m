w = randn(n0)/n0;  w = w-mean(w(:))+1/n0^2;
f = real(ifft2(fft2(f1).*repmat(fft2(w), [1 1 size(f0,3)])));
clf;
u = ones(n0,n0,3);
u(end/2-n/2+1:end/2+n/2, end/2-n/2+1:end/2+n/2, :) = f0; 
u(1,1,:)=0; u(2,1,:)=1;
imageplot(clamp(u), 'Input', 1,2,1);
u = f; u(1,1,:)=0; u(2,1,:)=1;
imageplot(clamp(u), 'Synthesized', 1,2,2);
