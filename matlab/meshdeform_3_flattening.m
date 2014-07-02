%% Spectral Mesh Flattening 
% This tour explores 2D flattening of 3D surfaces using spectral methods.

perform_toolbox_installation('signal', 'general', 'graph');


%CMT
rep = 'results/mesh_flattening/';
if not(exist(rep))
    mkdir(rep);
end
%CMT

%% Spectral Mesh Flattening
% Mesh flattening finds 2D locations that minimize a variational energy
% with a non-degenracy constraint (for instance maximal variance). 

%% 
% For the Dirichlet energy (Sobolev norm), the
% resulting location are described by the first eigenvectors of the
% Laplacian. 

%%
% This method is refered to as "Laplacian eigenmaps" in manifold learning,
% see:

%%
%

%% 
% The advantage over fixed boundary harmonic parameterization is
% that the boundary of the flattened domain is not fixed, but the drawback
% is that the parameterization is not guaranteed to be valid (bijective).

%% 
% First load a mesh.

name = 'nefertiti'; 
options.name = name;
[vertex,faces] = read_mesh(name);
n = size(vertex,2);

%%
% Display it.

clf;
plot_mesh(vertex,faces, options);
shading faceted;

%CMT
saveas(gcf, [rep 'flattening-' name '-original.eps'], 'epsc');
%CMT

%%
% Compute the mesh Laplacian matrix.

options.symmetrize = 1;
options.normalize = 0;
L = compute_mesh_laplacian(vertex,faces,'conformal',options);

%% 
% Compute the eigenvalues and eigenvectors

[U,S] = eig(full(L)); S = diag(S);
[S,I] = sort(S,'ascend'); U = U(:,I);

%%
% The vertex positions are the eigenvectors 2 and 3.

vertexF = U(:,2:3)';

%%
% Use translation / rotation to align the parameterization.

icenter = 88;
irotate = 154;
vertexF = vertexF - repmat(vertexF(:,icenter), [1 n]);
theta = -pi/2+atan2(vertexF(2,irotate),vertexF(1,irotate));
vertexF = [vertexF(1,:)*cos(theta)+vertexF(2,:)*sin(theta); ...
           -vertexF(1,:)*sin(theta)+vertexF(2,:)*cos(theta)];

%%
% Display the flattened mesh. 

clf;
plot_mesh(vertexF,faces);

%CMT
saveas(gcf, [rep 'flattening-' name '-lapleig-conformal.eps'], 'epsc');
%CMT

%EXO
%% Perform the same flattening, but with the combinatorial Laplacian.
L = compute_mesh_laplacian(vertex,faces,'combinatorial',options);
[U,S] = eig(full(L)); S = diag(S);
[S,I] = sort(S,'ascend'); U = U(:,I);
vertexF = U(:,2:3)';
vertexF = vertexF - repmat(vertexF(:,icenter), [1 n]);
theta = -pi/2+atan2(vertexF(2,irotate),vertexF(1,irotate));
vertexF = [vertexF(1,:)*cos(theta)+vertexF(2,:)*sin(theta); ...
           -vertexF(1,:)*sin(theta)+vertexF(2,:)*cos(theta)];
clf;
plot_mesh(vertexF,faces);
%EXO

%CMT
saveas(gcf, [rep 'flattening-' name '-lapleig-combinatorial.eps'], 'epsc');
%CMT

%% Geodesic Embedding (Isomap)
% Another (nonlinear) embedding can be computed by minimizing the geodesic
% distortion between points on the surface and points over the
% parameterized domain. 

%%
% First we compute the geodesic distance on the mesh using the Fast Marching algorithm.

D = zeros(n);
for i=1:n
    D(:,i) = perform_fast_marching_mesh(vertex,faces,i);
end

%%
% Enforce symmetry.

D = (D+D')/2;

%%
% Compute the centered matrix.

J = eye(n) - ones(n)/n;
W = -J*(D.^2)*J;

%%
% Diagonalize the centered matrix.

[U,S] = eig(W);
S = diag(S);
[S,I] = sort(S,'descend'); U = U(:,I);

%%
% Display the decay of the eigenvalues. If the mesh was isometric to the
% plane, then only the two largest eigenvalues would be non zero.

clf;
hh = plot(S(1:30), '.-'); axis('tight');
set(hh, 'LineWidth', 2);

%%
% Isomap embedding is defined from the two largest eigenvalues.

vertexF = U(:,1:2)' .* repmat(sqrt(S(1:2)), [1 n]);

%%
% Align the parameters.

vertexF = vertexF - repmat(vertexF(:,icenter), [1 n]);
theta = -pi/2+atan2(vertexF(2,irotate),vertexF(1,irotate));
vertexF = [vertexF(1,:)*cos(theta)+vertexF(2,:)*sin(theta); ...
           -vertexF(1,:)*sin(theta)+vertexF(2,:)*cos(theta)];

%% 
% Display and compare with Laplacian embeddedding.

clf;
plot_mesh(vertexF,faces,options);


%CMT
saveas(gcf, [rep 'flattening-' name '-isomap-strain.eps'], 'epsc');
%CMT

%EXO
%% Compute the embedding using Stress minimization with SMACOF.
%% See the numerical tours on bending invariants for more details.
niter = 150;
stress = [];
vertexS = vertexF;
ndisp = [1 5 10 min(niter,100) Inf];
k = 1;
clf;
for i=1:niter
    if ndisp(k)==i
        subplot(2,2,k);
        plot_mesh(vertexS,faces,options);
        k = k+1;
    end
    % Compute the distance matrix.
    D1 = repmat(sum(vertexS.^2,1),n,1);
    D1 = sqrt(D1 + D1' - 2*vertexS'*vertexS);
    % Compute the scaling matrix.
    B = -D./max(D1,1e-10);
    B = B - diag(sum(B));
    % update
    vertexS = (B*vertexS')' / n;
    % Xstress = Xstress-repmat(mean(Xstress,2), [1 n]);
    % record stress
    stress(end+1) = sqrt( sum( abs(D(:)-D1(:)).^2 ) / n^2 );
end
%EXO


%%
% Plot stress evolution during minimization.

clf;
plot(stress(2:end), '.-');
axis('tight');

%EXO
%% Compute mesh parameterization using a circle as boundary.
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
%EXO


%CMT
saveas(gcf, [rep 'flattening-' name '-parameterization.eps'], 'epsc');
%CMT
