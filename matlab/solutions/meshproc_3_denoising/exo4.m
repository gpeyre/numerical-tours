ntests = 15;
muList = linspace(3, 15, ntests)/5;
errR = [];
for i=1:ntests
    mu = muList(i);
    A = speye(n,n)+mu*L;
    for k=1:3
        Xmu(k,:) = perform_cg(A,X(k,:)')';
    end
    errR(i) = snr(X0,Xmu);
end
clf;
plot(muList, errR, '.-'); axis('tight');
set_label('\mu', 'SNR');
