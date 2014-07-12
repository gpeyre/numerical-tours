sv = .2;
slist = linspace(1,8*2,4);
clf;
for i=1:4
    sx = slist(i);
    imageplot( bilateral_lin(f0,sx,sv), ['\sigma_x=' num2str(sx,2)], 2,2,i );
end
