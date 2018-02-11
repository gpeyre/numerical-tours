delta = zeros(N,N);
for i=1:N
%    progressbar(i,N);
    [delta(:,i),S,Q] = perform_fast_marching_mesh(V, F, i);
end
delta = (delta+delta)'/2;
