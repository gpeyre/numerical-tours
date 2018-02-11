niter = 5;
niter1 = 5;
M1 = randn(n,n,3);
M1 = perform_color_matching(M1,M,niter1);
% precompute the wavelet transform of the image
for c=1:3
    MW(:,:,:,c) = perform_wavelet_transf(M(:,:,c), Jmin, +1, options);
end
Msvg1 = {}; Msvg = {};
for k=1:niter
    % wavelet matching
    for c=1:3
        % transform
        MW1 = perform_wavelet_transf(M1(:,:,c), Jmin, +1, options);
        % equalize
        for i=1:size(MW1,3)
            MW1(:,:,i) = perform_hist_eq(MW1(:,:,i), MW(:,:,i,c));
        end
        M1(:,:,c) = perform_wavelet_transf(MW1, Jmin, -1, options);
    end
    % spacial matching
    M1 = perform_color_matching(M1,M,niter1);
    Msvg{end+1} = M1;
end
% Display.
clf;
imageplot(M, 'Exemplar', 1,2,1);
imageplot(M1, 'Synthesized', 1,2,2);
