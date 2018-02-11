mydisp = lambda x: np.log10(abs(pyl.fftshift(x,axes=(0,1))) + 1e-5)

plt.figure(figsize = (10,10))
imageplot(np.mean(mydisp(pyl.fft2(f, axes=(0,1))),2), 'Original', [1, 2, 1])
imageplot(np.mean(mydisp(pyl.fft2(p, axes=(0,1))),2), 'Periodic layer', [1, 2, 2])