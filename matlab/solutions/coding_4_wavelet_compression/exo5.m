Tlist = linspace(.03,.6,20);
err = []; nbits = [];
for i=1:length(Tlist)
    T = Tlist(i);
    % quantize
    MWI = floor(abs(MW/T)).*sign(MW);
    MWQ = sign(MWI) .* (abs(MWI)+.5) * T;
    % inverse
    MQ = perform_wavelet_transf(MWQ,Jmin, -1);
    % error
    err(i) = snr(M,MQ);
    % bits
    nbits(i) = compute_entropy(MWI(:));
end
clf;
hh = plot(nbits,err); axis('tight');
set_label('bpp','SNR');
if using_matlab()
    set(hh,'LineWidth',2);
end
