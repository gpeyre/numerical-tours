Tlist = linspace(.4,1.2,15)*sigma;
err = [];
for i=1:length(Tlist)
    ST = perform_thresholding(S, Tlist(i), 'soft');
    xT = perform_stft(ST,wlist,qlist, options);
    err(i) = snr(x0,xT);
end
clf;
plot(Tlist/sigma, err); axis('tight');
set_label('T/\sigma', 'SNR');
