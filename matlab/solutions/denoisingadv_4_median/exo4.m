klist = 1:4;
err = [];
for k=klist
    Mmed = perform_median_filtering(M,k);
    err(k) = snr(M0,Mmed);
end
clf;
plot(klist*2+1, err, '.-'); axis('tight');
set_label('w','SNR');
