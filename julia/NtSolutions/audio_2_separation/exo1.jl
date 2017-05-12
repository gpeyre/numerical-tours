for i in 1:p
    Y[:,:,i] = perform_stft(y[:,i],w,q,n)
    figure(figsize = (15,10))
    plot_spectrogram(Y[:,:,i],"Source #$i")
end
