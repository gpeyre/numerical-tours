% spacial matching
niter = 4;
M1 = randn(n);
clf;
for k=1:niter
    M1 = perform_hist_eq(M1, M);
    % wavelet matching
    MW1 = perform_wavelet_transf(M1, Jmin, +1, options);
    for i=1:size(MW,3)
        MW1(:,:,i) = perform_hist_eq(MW1(:,:,i), MW(:,:,i));
    end
    M1 = perform_wavelet_transf(MW1, Jmin, -1, options);
    % display
    imageplot(M1, strcat(['Iteration ' num2str(k)]), 2,2,k);
end
