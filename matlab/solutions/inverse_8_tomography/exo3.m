alist = linspace(0,.99*pi/2,4);
clf;
for i=1:length(alist)
    a = alist(i);
    imageplot(clamp(rotation(f0,a)), ['\theta=' num2str(a,2)], 2,2,i);
end
