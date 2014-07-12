options.h = compute_wavelet_filter('Daubechies',2);
MW = perform_wavortho_transf(M,0,+1,options);
MWT = perform_thresholding(MW(:), round(rho*n), 'strict');
MWT = reshape(MWT, size(M));
M1 = perform_wavortho_transf(MWT,0,-1,options);
clf;
imageplot(clamp(crop(M1)));
