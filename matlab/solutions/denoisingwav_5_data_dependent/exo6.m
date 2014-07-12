slist = [.1 .2 .3 .6];
clf;
for i=1:length(slist)
    Wu = gamrnd(1/slist(i)^2, slist(i)^2,n,n);
    imageplot(f0.*Wu, strcat(['\sigma=' num2str(slist(i))]), 2,2,i );
end
