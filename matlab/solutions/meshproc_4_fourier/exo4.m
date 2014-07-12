ntests = 50;
Tlist = linspace(.1,3,ntests);
err = [];
for i=1:ntests
    T = Tlist(i);
    % decoding
    pvertexI = floor(abs(pvertex/T)).*sign(pvertex);
    pvertexQ = sign(pvertexI) .* (abs(pvertexI)+.5) * T;
    vertex1 = pvertexQ*U';
    % entropic
    t = min(pvertexI(:)):max(pvertexI(:));
    h = hist( pvertexI(:), t );
    h = max(h,1e-10); h = h/sum(h);
    E = -sum( log2(h).*h );
    % recode
    nbits(i) = 3*E;
    err(i) = snr(vertex,vertex1);
end
% plot
clf;
plot(nbits, err); axis('tight');
set_label('nb.bits', 'SNR');
