ntests = 20;
Mlist = round( linspace(.02*n^2, .4*n^2,ntests)  );
su = [];
tic; x1 = K \ y;    t1 = toc();
ntrials = 14;
for i=1:ntests
    KWT = perform_thresholding(KW, Mlist(i), 'largest');
    KWT = sparse(KWT);
    T = 0;
    for k=1:ntrials
        tic; x1 = KWT \ yW; 
        T = T + toc();
    end
    T = T/ntrials;
    su(end+1) = 1-T/t1;
end
clf;
plot(Mlist/n^2, su, '.-'); axis('tight');
set_label('M/n^2', 'av. speed up');
