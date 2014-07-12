MW = perform_haar_transf(Mnoisy, 1, +1);
Tlist = linspace(1,4,20)*sigma;
err_hard = []; err_soft = [];
for i=1:length(Tlist)
    MWT = perform_thresholding(MW, Tlist(i), 'hard');
    M1 = perform_haar_transf(MWT, 1, -1);
    err_hard(i) = snr(M,M1);
    MWT = perform_thresholding(MW, Tlist(i), 'soft');
    M1 = perform_haar_transf(MWT, 1, -1);
    err_soft(i) = snr(M,M1);
    if i>1 & err_soft(i)>max(err_soft(1:i-1))
        Mwav = M1;
    end
end
clf;
plot(Tlist/sigma, [err_hard; err_soft]', '.-');
axis('tight');
set_label('T/sigma', 'SNR');
legend('hard', 'soft');
