Sn = perform_stft(xn,w,q, options);
Tlist = linspace(.8,2.5,20)*sigma;
err_hard = []; err_soft = [];
for i=1:length(Tlist);
    % soft thresholding
    SnT = perform_thresholding(Sn,Tlist(i)/2, 'soft');
    x1 = perform_stft(SnT, w,q, options);
    err_soft(i) = snr(x,x1);
    % hard thresholding
    SnT = perform_thresholding(Sn,Tlist(i), 'hard');
    x1 = perform_stft(SnT, w,q, options);
    err_hard(i) = snr(x,x1);
end
% retrieve best hard thresholding result
[snr_hard,i] = max(err_hard);
SnT = perform_thresholding(Sn,Tlist(i), 'hard');
x1 = perform_stft(SnT, w,q, options);
% display the error curves
clf;
plot(Tlist/sigma, [err_hard(:) err_soft(:)]);
axis('tight');
legend('Hard', 'Soft');
set_graphic_sizes([], 20,2);
set_label('T/sigma', 'SNR');
