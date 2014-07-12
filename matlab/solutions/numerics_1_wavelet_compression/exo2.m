ntests = 20;
Mlist = round( linspace(.02*n^2, .4*n^2,ntests)  );
err = [];
for i=1:ntests
    M = Mlist(i);
    KWT = perform_thresholding(KW, M, 'largest');
    xTW = KWT \ yW;
    xT = perform_wavelet_transf(xTW, Jmin, -1, options);
    err(end+1) = norm(x-xT) / norm(x);    
end
clf;
plot(Mlist/n^2, log10(err), '.-'); axis('tight');
set_label('M/n^2', 'log_{10}(|x-xT|/|x|)');
