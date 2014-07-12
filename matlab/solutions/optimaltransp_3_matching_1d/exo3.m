p = 100;
tlist = linspace(0,1,5);
for i=1:length(tlist)
    a = ft(tlist(i));
    subplot(length(tlist), 1, i);
    [h,t] = hist(a(:), p);
    bar(t,h*p/n^2);
    axis([0 1 0 6]);
end
%EXO