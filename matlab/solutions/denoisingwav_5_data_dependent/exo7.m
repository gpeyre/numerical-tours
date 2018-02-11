ntest = 30;
sigma = std(f(:)-f0(:));
Tlist = linspace(2.8,4.5,ntest)*sigma;
fW = perform_wavelet_transf(f,Jmin,+1,options);
err = [];
for i=1:ntest
    fWT = perform_thresholding(fW, Tlist(i), 'hard');
    f1 = perform_wavelet_transf(fWT,Jmin,-1,options);    
    err(i) = snr(f0,f1);
end
% display
clf;
h = plot(Tlist/sigma,err);
axis('tight');
if using_matlab()
    set(h, 'LineWidth', 2);
end
set_label('T/sigma', 'SNR');
% best result
[tmp,i] = max(err);
fWT = perform_thresholding(fW, Tlist(i), 'hard');
fMult = perform_wavelet_transf(fWT,Jmin,-1,options);
