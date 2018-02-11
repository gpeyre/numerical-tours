% compute laplacian
W0 = sparse(n,n);
for i=1:3
   i1 = mod(i-1,3)+1;
   i2 = mod(i  ,3)+1;
   i3 = mod(i+1,3)+1;
   pp = vertex0(:,faces(i2,:)) - vertex0(:,faces(i1,:));
   qq = vertex0(:,faces(i3,:)) - vertex0(:,faces(i1,:));
   pp = pp ./ repmat( sqrt(sum(pp.^2,1)), [3 1] );
   qq = qq ./ repmat( sqrt(sum(qq.^2,1)), [3 1] );
   ang = acos(sum(pp.*qq,1));
   u = cot(ang);
   u = clamp(u, 0.01,100);
   W0 = W0 + sparse(faces(i2,:),faces(i3,:),u,n,n);
   W0 = W0 + sparse(faces(i3,:),faces(i2,:),u,n,n);
end
d = full( sum(W0,1) );
D = spdiags(d(:), 0, n,n);
L0 = D - W0;
% 
Delta0 = zeros(3,n);
Delta0(3,I1) = 1;
% laplacian 
L01 = L0; L01(I,:) = 0; L01(I + (I-1)*n) = 1;
% deform
vertex1 = vertex0 + ( L01 \ Delta0' )';
% Display it.
clf;
plot_mesh(vertex1,faces);
view(-150,45);
zoom(.8);
