slist = [4 6 10 15];
clf;
for i=1:length(slist)
    sigma = slist(i);
    subplot( 2,2,i );
    plot_levelset( delta(blur(f0,sigma)) ,0,f0);
    title(['\sigma=' num2str(sigma)]);
end
