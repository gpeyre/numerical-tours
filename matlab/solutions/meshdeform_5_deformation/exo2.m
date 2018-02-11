J = reshape(1:p^2,p,p);
q = round(.3*p); sel = q:p-q+1;
I0 = unique( [J(1,:) J(end,:) J(:,1)' J(:,end)']' );
I1 = unique( [J(q,sel) J(p-q+1,sel) J(sel,q)' J(sel,p-q+1)']' );
I = [I0; I1];
% 
Delta0 = zeros(3,n);
Delta0(3,I1) = 1;
% laplacian 
L1 = L; L1(I,:) = 0; L1(I + (I-1)*n) = 1;
% deform
vertex1 = vertex+( L1 \ Delta0' )';
% Display it.
clf;
plot_mesh(vertex1,faces);
view(-150,45);
zoom(.8);
