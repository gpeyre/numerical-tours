W = make_sparse(n,n);
for i=1:3
   i1 = mod(i-1,3)+1;
   i2 = mod(i  ,3)+1;
   i3 = mod(i+1,3)+1;
   pp = X(:,F(i2,:)) - X(:,F(i1,:));
   qq = X(:,F(i3,:)) - X(:,F(i1,:));
   % normalize the vectors   
   pp = pp ./ repmat( sqrt(sum(pp.^2,1)), [3 1] );
   qq = qq ./ repmat( sqrt(sum(qq.^2,1)), [3 1] );
   % compute angles
   a = 1 ./ tan( acos(sum(pp.*qq,1)) );
   a = max(a, 1e-2); % avoid degeneracy
   W = W + make_sparse(F(i2,:),F(i3,:), a, n, n );
   W = W + make_sparse(F(i3,:),F(i2,:), a, n, n );
end
% Compute the symmetric Laplacian matrix.
d = full( sum(W,1) );
D = spdiags(d(:), 0, n,n);
L = D - W;
%
options.verb = 0;
B = compute_boundary(F, options);
%
L1 = L;
L1(B,:) = 0;
for i=1:length(B)
    L1(B(i),B(i)) = 1;
end
