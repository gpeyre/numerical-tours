k = round(.8*n); k = round(k/2)*2; % even number
Mf = fft2(M);
Mf(n/2-k/2+2:n/2+k/2, n/2-k/2+2:n/2+k/2) = 0;
Mh = real( ifft2(Mf) );
% display
clf;
imageplot( crop(M), 'Image', 1,2,1);
imageplot(clamp( crop(Mh)), 'Low pass filtered', 1,2,2);
