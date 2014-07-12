err = [];
kmax = 12;
for k=1:kmax
    fmed = x;
    for x = 1:n
        sel = x-k:x+k;
        sel = mod(sel-1,n)+1;
        fmed(x) = median(f(sel));
    end
    err(k) = snr(f0,fmed);
end
clf;
plot(1:kmax, err, '.-'); axis('tight');
set_label('w', 'SNR');
