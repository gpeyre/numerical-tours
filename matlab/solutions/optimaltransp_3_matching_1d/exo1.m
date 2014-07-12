plist = [10 30 100 1000];
clf;
for i=1:length(plist)
    subplot(length(plist), 1, i);
    [h,t] = hist(f(:), plist(i));
    bar(t,h*plist(i)/n^2);
    axis([0 1 0 2.5]);
end
