function plot_spectrogram(S,f, options)

// plot_spectrogram - display the spectrogram
//
//   plot_spectrogram(S,f, options);
//
//   f is the 1D signal.
//   S is the spectrogram.
//
//  Copyright (c) 2008 Gabriel Peyre

// keep only low frequencies.
S = rescale(abs(S(1:$/2,:)));

options.null = 0;
if argn(2)<2
    f = [];
end
normalization = getoptions(options , 'normalization', 'log');
if strcmp(normalization, 'log')
    S = log(S+1e-4);
end

if not(isempty(f))
    n = size(f,1);
    clf;
    subplot(2,1,1);
    plot(1:n,f); axis tight;
    subplot(2,1,2);
end
imageplot(S);

set_colormap('jet');

//axis normal; 
//axis on;
//colormap jet(256);

endfunction