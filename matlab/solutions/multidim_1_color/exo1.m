rgb2yuv = @(f)applymat(f,T^(-1));
f1 = rgb2yuv(U1);
clf;
imageplot(f, 'Image', 1,2,1);
imageplot(clamp(f1), 'Modified', 1,2,2);
