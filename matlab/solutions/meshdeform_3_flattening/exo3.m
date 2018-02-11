% Cotan weights
options.symmetrize = 1;
options.normalize = 0;
L = compute_mesh_laplacian(vertex,faces,'conformal',options);
% boundary
options.verb = 0;
boundary = compute_boundary(faces, options);
% fixed positions
p = length(boundary);
t = linspace(0,2*pi,p+1)'; t(p) = [];
x0 = cos(t); y0 = sin(t);
% system
L1 = L;
L1(boundary,:) = 0;
L1(boundary + (boundary-1)*n) = 1;
% Set up the right hand sizes with the fixed position.
Rx = zeros(n,1); Rx(boundary) = x0;
Ry = zeros(n,1); Ry(boundary) = y0;
% solve
x = L1 \ Rx;
y = L1 \ Ry;
vertexF = [x';y'];
% align
vertexF = vertexF - repmat(vertexF(:,icenter), [1 n]);
theta = -pi/2+atan2(vertexF(2,irotate),vertexF(1,irotate));
vertexF = [vertexF(1,:)*cos(theta)+vertexF(2,:)*sin(theta); ...
           -vertexF(1,:)*sin(theta)+vertexF(2,:)*cos(theta)];
% display
clf;
plot_mesh(vertexF,faces,options);
