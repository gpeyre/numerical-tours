fwLow = zeros(n);
fwLow(1:2^J,1:2^J) = fw(1:2^J,1:2^J);
fLow = WavI(fwLow);
myplot = @(f1)imageplot(clamp(f1), ['PSNR=' num2str(psnr(f,f1),3) 'dB']);
clf;
myplot(fLow);
