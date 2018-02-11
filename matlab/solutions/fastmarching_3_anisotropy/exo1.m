npoints = 10;
nstart = floor( rand(2,npoints)*n ) + 1;
aniso_list = [.02 .1 .2 .5];
clf;
for i=1:length(aniso_list)
    anisotropy = aniso_list(i);
    H = perform_tensor_recomp(e1,e2, ones(n),ones(n)*1/anisotropy );
    [D, dUx, dUy, Vor, L] = fm2dAniso([hx;hy], H, pstart);
    imageplot(convert_distance_color(D), ['Anisotropy=' num2str(anisotropy)],2,2,i);
end
