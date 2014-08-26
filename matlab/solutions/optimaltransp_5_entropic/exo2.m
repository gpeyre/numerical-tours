glist = [.1 .01 .001 .0001];
niter = 500;
clf;
for ig=1:length(glist)
    gamma = glist(ig);
    pi = exp( -C/gamma );
    for i=1:niter
        pi = ProjC2( ProjC1(pi,p), q);
    end
    imageplot(normalizeMax(pi), ['\gamma=' num2str(gamma)], 2,2,ig);
end
