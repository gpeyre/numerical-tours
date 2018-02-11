alist = linspace(0,1,4);
clf;
for i=1:length(alist)
    a = alist(i);
    imageplot(clamp(shearx(f0,a)), ['a=' num2str(a,2)], 2,2,i);
end
