% display
clf;
subplot(2,1,1);
plot(xn); axis([1 n -1.2 1.2]);
set_graphic_sizes([], 20);
title(strcat(['Noisy signal, SNR=', num2str(snr(x,xn),4), 'dB']));
subplot(2,1,2);
plot(x1); axis([1 n -1.2 1.2]);
set_graphic_sizes([], 20);
title(strcat(['Denoised signal, SNR=', num2str(snr_hard,4), 'dB']));
% hear
sound(x1,fs);
