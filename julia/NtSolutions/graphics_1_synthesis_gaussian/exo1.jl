mydisp = x-> log10(abs(fftshift(x) + 1e-5))

figure(figsize = (10,10))
imageplot(mean(mydisp(fft(f,(1,2))),3)[:,:,1], "Original", [1, 2, 1])
imageplot(mean(mydisp(fft(p,(1,2))),3)[:,:,1], "Periodic layer", [1, 2, 2]);
