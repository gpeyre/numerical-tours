err = [];
shift_list = 1:7;
T = 3*sigma;
for m=shift_list
    [dY,dX] = meshgrid(0:m-1,0:m-1);
    fTI = zeros(n,n);
    for i=1:m^2
        fS = circshift(f,[dX(i) dY(i)]);
        a = perform_wavelet_transf(fS,Jmin,1,options);
        aT = perform_thresholding(a,T,'hard');
        fS = perform_wavelet_transf(aT,Jmin,-1,options);
        fS = circshift(fS,-[dX(i) dY(i)]);
        fTI = (i-1)/i*fTI + 1/i*fS;
    end
    err(m) = snr(f0,fTI);
end
clf;
h = plot(shift_list, err, '.-');
if using_matlab()
    set(h, 'LineWidth', 2);
end
axis('tight');
set_label('m', 'SNR');
