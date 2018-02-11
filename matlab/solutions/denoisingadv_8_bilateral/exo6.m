sx = 4*2;
slist = linspace(.05,.4,4);
clf;
for i=1:4
    sv = slist(i);
    imageplot( bilateral_lin(f0,sx,sv), ['\sigma_v=' num2str(sv,2)], 2,2,i );
end
