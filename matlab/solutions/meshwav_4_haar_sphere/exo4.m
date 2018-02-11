% Forward transform
rho = .1;
rho = .05;
nj = length(f);
fw = f;
for j=1:J-1
    fj = fw(1:nj);
    fj = reshape(fj, [nj/4 4]);
    fj = fj*U;
    fw(1:nj) = fj(:);
    nj = nj/4;
end
% thresholding
fwT = perform_thresholding(fw, round(rho*n), 'strict');
% Backward transform
nj = size(face{1},2);
f1 = fwT;
for j=1:J-1
    fj = f1(1:4*nj);
    fj = reshape(fj, [nj 4]);
    fj = fj*U';
    f1(1:4*nj) = fj(:);
    nj = nj*4;
end
% display
clf;
options.face_vertex_color = clamp(f1);
plot_mesh(vertex{end}, face{end}, options);
view(vv);
colormap gray(256);
lighting none;
