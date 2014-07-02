function plot_spectrogram(S,f, options)

% plot_spectrogram - display the spectrogram
%
%   plot_spectrogram(S,f, options);
%
%   f is the 1D signal.
%   S is the spectrogram.

options.null = 0;
if nargin<2
    f = [];
end

if iscell(S)
    e = 0;
    if ~isempty(f)
        e=1;
    end
    clf;
    if e==1
        subplot(length(S)+e,1,1);
        plot(1:length(f),f); axis tight;        
    end
    for i=1:length(S)
        subplot(length(S)+e,1,i+e);
        plot_spectrogram(S{i}, [], options);
    end
    return;
end

% keep only low frequencies.
S = rescale(abs(S(1:end/2,:)));

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
axis normal; 
axis on;
colormap jet(256);