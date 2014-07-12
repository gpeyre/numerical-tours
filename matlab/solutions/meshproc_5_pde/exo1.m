W = make_sparse(n,n);
for i=1:3
   i1 = mod(i-1,3)+1;
   i2 = mod(i  ,3)+1;
   i3 = mod(i+1,3)+1;
   pp = X0(:,F(i2,:)) - X0(:,F(i1,:));
   qq = X0(:,F(i3,:)) - X0(:,F(i1,:));
   % normalize the vectors   
   pp = pp ./ repmat( sqrt(sum(pp.^2,1)), [3 1] );
   qq = qq ./ repmat( sqrt(sum(qq.^2,1)), [3 1] );
   % compute angles
   ang = acos(sum(pp.*qq,1));
   W = W + make_sparse(F(i2,:),F(i3,:), 1 ./ tan(ang), n, n );
   W = W + make_sparse(F(i3,:),F(i2,:), 1 ./ tan(ang), n, n );
end
