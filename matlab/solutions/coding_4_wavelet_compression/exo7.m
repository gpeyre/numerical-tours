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
disp(['nb.bis, whole:    ' num2str(Ewhole,3) ' bpp']);
disp(['nb.bis, separate: ' num2str(Esep/n^2,3) ' bpp']);
