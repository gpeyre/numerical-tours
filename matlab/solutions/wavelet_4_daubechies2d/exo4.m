% number of kept coefficients
m = round(n^2/16); 
% compute the threshold T
Jmin = 1;
fW = perform_wavortho_transf(f,Jmin,+1); 
% select threshold
v = sort(abs(fW(:))); 
if v(1)<v(n^2)
    v = reverse(v);
end
% inverse transform
T = v(m);
fWT = fW .* (abs(fW)>=T);
% inverse
Mnlin = perform_wavortho_transf(fWT,Jmin,-1);
% display
clf;
imageplot(clamp(fLin),strcat(['Linear, SNR=' num2str(snr(f,fLin),3) 'dB']), 1,2,1 );
imageplot(clamp(Mnlin),strcat(['Non-linear, SNR=' num2str(snr(f,Mnlin),3) 'dB']), 1,2,2 );
