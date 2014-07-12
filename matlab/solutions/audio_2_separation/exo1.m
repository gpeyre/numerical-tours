clf;
for i=1:p
    Y(:,:,i) = perform_stft(y(:,i),w,q, options);
    subplot(p,1,i);
    plot_spectrogram(Y(:,:,i));
    set_graphic_sizes([], 20);
    title(strcat('Micro 1',num2str(i)));
end
