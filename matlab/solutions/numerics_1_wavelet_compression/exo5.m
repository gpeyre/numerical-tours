Mlist = [.01 .05 .1]*n^2 / 20;
ntrials = 20;
clf;
for i=1:length(Mlist)
    M = Mlist(i);
    KWT = perform_thresholding(KW, M, 'largest');
    KWT = sparse(KWT);
    t2 = 0;
    for k=1:ntrials
        tic;
        xW = perform_wavelet_transf(x, Jmin, +1, options);
        yWT = KWT * xW;
        yT = perform_wavelet_transf(yWT, Jmin, -1, options);
        t2 = t2 + toc();
    end
    t2 = t2 / ntrials;
    subplot(length(Mlist), 1, i);
    plot(yT); axis('tight');
    title(['M/n^2=' num2str(M/n^2) ', SNR=' num2str(snr(y,yT), 3) ', speed up=' num2str(1-t2/t0, 3) '.']);
end
%EXO