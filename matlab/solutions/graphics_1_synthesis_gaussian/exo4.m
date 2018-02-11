clf;
for i=1:4
    w = randn(n)/n;  w = w-mean(w(:))+1/n^2;
    f = real(ifft2(fft2(f0).*repmat(fft2(w), [1 1 size(f0,3)])));
    u = f; u(1,1,:)=0; u(2,1,:)=1;
    subplot(2,2,i);
    imageplot(clamp(u));
end
