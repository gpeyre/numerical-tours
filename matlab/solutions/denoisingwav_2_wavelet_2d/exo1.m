Tlist = linspace(.8,4.5,25)*sigma;
err_soft = []; err_hard = [];
for i=1:length(Tlist)
    aT = perform_thresholding(a,Tlist(i),'hard');
	fWav = perform_wavelet_transf(aT,Jmin,-1,options);
    err_hard(i) = snr(f0,fWav);
    aT = perform_thresholding(a,Tlist(i),'soft');
    aT(1:2^Jmin,1:2^Jmin) = a(1:2^Jmin,1:2^Jmin);
	fWav = perform_wavelet_transf(aT,Jmin,-1,options);
    err_soft(i) = snr(f0,fWav);
end
clf;
h = plot(Tlist/sigma, [err_hard(:) err_soft(:)]); axis('tight');
if using_matlab()
    set(h, 'LineWidth', 2);
end
set_label('T/\sigma', 'SNR');
legend('Hard', 'Soft');
