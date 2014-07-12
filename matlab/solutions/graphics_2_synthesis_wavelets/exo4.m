M1 = randn(n,n,3);
niter1 = 200;
delta = [];
for k=1:niter1
    [U,R] = qr(randn(3));
    d = reshape(M,[n^2 3])*U;
    d1 = reshape(M1,[n^2 3])*U;
    delta(k) = 0;
    for c=1:3
        [tmp,I] = sort(d(:,c)); 
        [tmp,I1] = sort(d1(:,c)); 
        delta(k) = delta(k) + norm(d1(I1,c) - d(I,c))^2;
        d1(I1,c) = d(I,c);
        % d1(:,c) = perform_hist_eq(d1(:,c),d(:,c));
    end
    M1old = M1;
    M1 = reshape(d1*U',[n n 3]);
    % delta(k) = norm(M1(:)-M1old(:));
end
clf;
plot(log10(delta/norm(M1(:))), '.-');
axis('tight');
