lmin = 1;
lmax = [5 10 50 100];
clf;
for i=1:length(lmax)
    f1 = poissrnd( floor( rescale(f0u,lmin,lmax(i)) ) );
    imageplot(f1, strcat(['\lambda_{max}=' num2str(lmax(i))]), 2,2,i );
end
