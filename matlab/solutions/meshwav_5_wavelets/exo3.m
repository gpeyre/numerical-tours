clf;
nverts = size(vertex{end}, 2);
i = 0;
for j = [J-3 J-2]
    i = i+1;
    nj = size(vertex{j},2); nj1 = size(vertex{j+1},2);
    sel = nj+1:nj1-1;
    d = sum( abs(vertex{end}(:,sel)) );
    [tmp,k] = min(d); k = sel(k);
    fw2 = zeros(nverts,1); fw2(k) = 1;
    f2 = perform_wavelet_mesh_transform(vertex,face, fw2, -1, options);
    options.color = 'wavelets';
    options.use_color = 1;
    options.rho = .6;
    options.use_elevation = 1;
    options.view_param = [104,-40];
    subplot(1,2,i);
    plot_spherical_function(-vertex{end},face{end},f2, options); axis tight;
    title(['Wavelet, scale ' num2str(j)]);
end
options.view_param = [];
