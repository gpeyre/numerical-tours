T = 3*sigma;
w = 4;
[dX,dY,dZ] = ndgrid(0:w-1,0:w-1,0:w-1);
Mspin = zeros(n,n,n);
for i=1:w^3
    MnoisyC = circshift(Mnoisy, [dX(i) dY(i) dZ(i)]);
    % denoise
    MW = perform_wavortho_transf(MnoisyC, 1, +1, options);
    MWT = perform_thresholding(MW, T, 'hard');
    M1 = perform_wavortho_transf(MWT, 1, -1, options);
    % back
    M1 = circshift(M1, -[dX(i) dY(i) dZ(i)]);
    Mspin = Mspin*(i-1)/i + M1/i;
end
