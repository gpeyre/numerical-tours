fNN = bilateral_nn(f0,sx,sv);
fLin = bilateral_lin(f0,sx,sv);
c = [.5 .4]*n; q = 200;
clf;
imageplot(crop(fNN,q,c), 'Nearest neighbor', 1,2,1);
imageplot(crop(fLin,q,c), 'Linear', 1,2,2);
