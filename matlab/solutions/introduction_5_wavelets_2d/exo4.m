q = 30;
Tlist = linspace(2.5,4,q)*sigma;
err = [];
for i=1:q
    fWT = Thresh(fW,Tlist(i));    
    f1 = PsiS(fWT);
    err(i) = snr(f,f1);
end
clf;
plot(Tlist/sigma, err);
axis('tight');
set_label('T/\sigma', 'SNR');
