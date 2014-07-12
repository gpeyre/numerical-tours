ntest = 30;
sigma = 1;
Tlist = linspace(2,3.5,ntest)*sigma;
fW = perform_wavelet_transf( 2*sqrt(f+3/8) ,Jmin,+1,options);
err = [];
for i=1:ntest
    fWT = perform_thresholding(fW, Tlist(i), 'hard');
    f1 = perform_wavelet_transf(fWT,Jmin,-1,options);    
    % undo VST
    f1 = (f1/2).^2 - 3/8;
    % record error
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
f1 = perform_wavelet_transf(fWT,Jmin,-1,options);
fVST = (f1/2).^2 - 3/8;
