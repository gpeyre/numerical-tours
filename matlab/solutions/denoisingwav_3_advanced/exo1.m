mulist = linspace(1,12, 17);
nmu = length(mulist);
Tlist = linspace(.5,3.5,10)*sigma;
err = [];
for i=1:length(mulist)
    for j=1:length(Tlist)
        T = Tlist(j);
        MWT = perform_thresholding(MW,[T mulist(i)*T],'semisoft');
        MWT(1:2^Jmin,1:2^Jmin) = MW(1:2^Jmin,1:2^Jmin);
        MT = perform_wavelet_transf(MWT,Jmin,-1);
        err(i,j) = snr(M0,MT);
    end
end
clf;
imageplot(err); 
% Tlist/sigma, mulist, 
axis('tight');
set_label('T/\sigma', '\mu');
title('SNR of semi-soft thresholding');
axis('on');
