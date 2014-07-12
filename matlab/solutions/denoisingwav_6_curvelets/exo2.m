T = Tlist(i);
s = 4;
[dY,dX] = meshgrid(0:s-1,0:s-1);
Mcurv = zeros(n);
for i=1:s^2
    Ms = circshift(M, [dX(i) dY(i)]);
    MW = perform_curvelet_transform(Ms, options);
    MWT = perform_thresholding(MW, T, 'hard');
    Ms = perform_curvelet_transform(MWT, options);
    Ms = circshift(Ms, -[dX(i) dY(i)]);
    Mcurv = (1-1/i)*Mcurv + 1/i*Ms;
end
