fmed = x;
for x = 1:n
    sel = x-k:x+k;
    sel = mod(sel-1,n)+1;
    fmed(x) = median(f(sel));
end
clf;
subplot(2,1,1);
plot(clamp(f)); axis('tight');
title('Noisy signal');
subplot(2,1,2);
plot(clamp(fmed)); axis('tight');
title('Denoised signal');
