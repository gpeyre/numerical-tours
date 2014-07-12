tlist = linspace(0.5, 10, 4);
clf;
for i=1:length(tlist)
    t = tlist(i);
    imageplot(heat(f0,t), ['t=' num2str(t,2)], 2,2,i);
end
