mydisp = @(x)log10(abs(fftshift(x))+1e-5);
clf;
imageplot(mydisp(fft2(f)), 'Original',1,2,1);
imageplot(mydisp(fft2(p)), 'Periodic layer',1,2,2);
