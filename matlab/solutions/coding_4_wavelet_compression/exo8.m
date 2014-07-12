nbits1 = [];
for i=1:length(Tlist)
    T = Tlist(i);
    % quantize
    MWI = floor(abs(MW/T)).*sign(MW);
    % bits    
    Esep = 0;
    Jmax = log2(n)-1; Jmin = 4;
    for j = Jmax:-1:Jmin
        for q=1:3
            [selx,sely] = compute_quadsel(j,q);
            MWj = MWI(selx,sely);
            Esep = Esep + prod(size(MWj))*compute_entropy(MWj);
        end
    end
    Esep = Esep + prod(size(MWj))*compute_entropy(MWI(1:2^j,1:2^j));
    Ewhole = compute_entropy(MWI);
    nbits1(i) = Esep/n^2;
end
clf;
hh = plot([nbits(:)';nbits1(:)']',[err(:)';err(:)']'); axis('tight');
set_label('bpp','SNR');
legend('Whole','Separate');
if using_matlab()
    set(hh,'LineWidth',2);
end
