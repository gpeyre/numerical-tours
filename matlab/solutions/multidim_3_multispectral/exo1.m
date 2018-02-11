M1 = UWT;
for i=1:p
    M1(:,:,i) = perform_wavelet_transf(M1(:,:,i), Jmin, -1);
end
M1 = reshape( M1, [n*n p] )';
M1 = idct(M1);
M1 = reshape(M1', [n n p]);
%approximation error 
e = snr(M,M1);
% display
clf;
imageplot(M(:,:,rgbsel), 'Original', 1,2,1);
imageplot(clamp(M1(:,:,rgbsel)), strcat(['Approximated, SNR=' num2str(e,2) 'dB']), 1,2,2);
