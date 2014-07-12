alphalist = [.1 .3 .6 1];
for i=1:4
   alpha = alphalist(i);
   F = real( ifft2( fft2(f)./S.^alpha ) );
   subplot(2,2,i);
   imageplot(F, ['alpha=' num2str(alpha)]);
end
