# number of kept coefficients
m = round(n^2/16); 
# compute the threshold T
Jmin = 1;
fW = perform_wavortho_transf(f,Jmin,+1, h); 
# select threshold
v = sort(abs(fW[:])); 
if v[1]<v[n^2]
    v = reverse(v);
end
# inverse transform
T = v[m];
fWT = fW .* (abs(fW).>T);
# inverse
fnlin = perform_wavortho_transf(fWT,Jmin,-1, h);
# display
clf;
u1 = @sprintf("Linear, SNR=%.3fdB", snr(f,fLin));
u2 = @sprintf("Non-linear, SNR=%.3fdB", snr(f,fnlin)); 
imageplot(clamp(fLin),u1, 1,2,1 );
imageplot(clamp(fnlin),u2, 1,2,2 );
