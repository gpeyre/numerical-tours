E = zeros(n,n,nb_samples);
for i=1:nb_samples
    % progressbar(i,nb_samples);
    [d,tmp] = perform_fast_marching(W, samples(:,i), options);
    d(M==0) = 0;
    d(d==Inf) = 0;
    d(d>1e5) = 0;
    E(:,:,i) = d;
end
