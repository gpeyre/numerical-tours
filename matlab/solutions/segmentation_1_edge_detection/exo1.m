slist = [1 2 5 10];
clf;
for i=1:length(slist)
    sigma = slist(i);
    imageplot( blur(f0,sigma), ['\sigma=' num2str(sigma)], 2,2,i );
end
