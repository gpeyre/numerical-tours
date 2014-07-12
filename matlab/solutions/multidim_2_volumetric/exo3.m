ntests = 20;
slist = linspace(.01,1.5,ntests);
err = [];
for i=1:ntests
    h = exp( -(X.^2 + Y.^2 + Z.^2)/(2*slist(i)^2) );
    h = h/sum(h(:));
    Mh = real( ifftn( fftn(Mnoisy) .* fftn(fftshift(h)) ) );
    err(i) = snr(M,Mh);
    if i>1 && err(i)>max(err(1:i-1))
        Mblur = Mh;
    end
end
clf;
plot(slist,err, '.-');
axis('tight');
set_label('s', 'SNR');
