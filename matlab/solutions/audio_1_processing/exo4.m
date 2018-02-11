Sn = perform_stft(xn,w,q, options);
Tlist = linspace(1,3,6)/2*sigma;
snr_block = [];
bslist = 1:4;
[bsY,bsX] = meshgrid(bslist,bslist);
for k=1:length(bsX(:))
%    progressbar(k,length(bsX(:)));
    options.block_size = [bsX(k) bsY(k)];
    err = [];
    for i=1:length(Tlist);
        % hard thresholding
        SnT = perform_thresholding(Sn,Tlist(i), 'block',options);
        x1 = perform_stft(SnT, w,q, options);
        err(i) = snr(x,x1);
    end
    % retrieve best hard thresholding result
    [snr_block(k),t] = max(err);
    if t==1 | t==length(Tlist)
        warning('Out of bound reached');
    end
    Topt(k) = Tlist(t);
end
snr_block = reshape(snr_block, size(bsX));
clf;
imageplot(snr_block, 'SNR for several block sizes');
set_label('X block size','Y block size');
