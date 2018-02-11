clf;
slist = [.5 3 6 10];
for i=1:4
    s = slist(i);
    imageplot( Filter(f0,s), ['\sigma=' num2str(s)], 2,2,i );
end
