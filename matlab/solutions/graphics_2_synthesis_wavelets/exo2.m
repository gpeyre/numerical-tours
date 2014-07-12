niter = 3;
M1 = randn(n,n,3);
for c=1:3
    MW = perform_wavelet_transf(M(:,:,c), Jmin, +1, options);
    for k=1:niter
        M1(:,:,c) = perform_hist_eq(M1(:,:,c), M(:,:,c));
        % wavelet matching
        MW1 = perform_wavelet_transf(M1(:,:,c), Jmin, +1, options);
        for i=1:size(MW,3)
            MW1(:,:,i) = perform_hist_eq(MW1(:,:,i), MW(:,:,i));
        end
        M1(:,:,c) = perform_wavelet_transf(MW1, Jmin, -1, options);
    end
end
% Display.
clf;
imageplot(M, 'Exemplar', 1,2,1);
imageplot(M1, 'Synthesized', 1,2,2);
