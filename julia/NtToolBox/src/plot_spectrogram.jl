function plot_spectrogram(S,titre)
    S = abs(S[1:Int(size(S)[1]/2),:])
    S = log.(S + 1e-4)
    imshow(S, cmap = get_cmap("jet"), interpolation = "nearest")
    title(titre)
end
