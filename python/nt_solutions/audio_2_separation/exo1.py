for i in range(p):
    Y[:,:,i] = perform_stft(y[:,i],w,q,n)
    plt.figure(figsize = (15,10))
    plot_spectrogram(Y[:,:,i],"Source #%i" %(i+1))