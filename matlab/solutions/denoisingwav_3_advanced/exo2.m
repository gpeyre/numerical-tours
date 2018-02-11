Tlist = linspace(.5,4,20)*sigma;
err_hard = [];
err_soft = [];
err_stein = [];
for j=1:length(Tlist)
    T = Tlist(j);
    MWT = perform_thresholding(MW,T,'hard');
    MWT(1:2^Jmin,1:2^Jmin) = MW(1:2^Jmin,1:2^Jmin);
    MT = perform_wavelet_transf(MWT,Jmin,-1);
    err_hard(j) = snr(M0,MT);
    MWT = perform_thresholding(MW,T,'soft');
    MWT(1:2^Jmin,1:2^Jmin) = MW(1:2^Jmin,1:2^Jmin);
    MT = perform_wavelet_transf(MWT,Jmin,-1);
    err_soft(j) = snr(M0,MT);
    MWT = perform_thresholding(MW,T,'stein');
    MWT(1:2^Jmin,1:2^Jmin) = MW(1:2^Jmin,1:2^Jmin);
    MT = perform_wavelet_transf(MWT,Jmin,-1);
    err_stein(j) = snr(M0,MT);
end
clf;
hold('on');
plot(Tlist/sigma, err_hard, 'r');
plot(Tlist/sigma, err_soft, 'b');
plot(Tlist/sigma, err_stein, 'g');
hold('off');
axis('tight');
set_label('T/\sigma', 'SNR');
legend('Hard', 'Soft', 'Stein');
axis('on');
