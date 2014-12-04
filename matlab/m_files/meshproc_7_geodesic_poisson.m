%% Geodesic Distance with Poisson Equation
% This tour explores the method detailed in the paper <#biblio [CraneWeischedelWardetzky13]> to compute geodesic distance by 
% simply solving a Poisson equation.

perform_toolbox_installation('signal', 'general', 'graph');

%% Gradient, Divergence and Laplacian on Surfaces
% The topology of a triangulation is defined via a set of indexes \(\Vv = \{1,\ldots,n\}\)
% that indexes the \(n\) vertices, a set of edges \(\Ee \subset \Vv \times \Vv\)
% and a set of \(m\) faces \(\Ff \subset \Vv  \times \Vv \times \Vv\).

%%
% We load a mesh. The set of faces \(\Ff\) is stored in a matrix \(F \in
% \{1,\ldots,n\}^{3 \times m}\).
% The positions \(x_i \in \RR^3\), for \(i \in V\), of the \(n\) vertices
% are stored in a matrix \(X_0 = (x_{0,i})_{i=1}^n \in \RR^{3 \times n}\).

clear options;
name = 'elephant-50kv';
options.name = name; % useful for displaying
[X,F] = read_mesh(name);

%% 
% Number \(n\) of vertices and number \(m\) of faces.

n = size(X,2);
m = size(F,2);

%%
% Display the mesh in 3-D.

options.lighting = 1;
clf;
plot_mesh(X,F,options); 
axis('tight');

%%
% In this tour, we use piecewise linear finite element to compute the gradient
% operators, which in turns allows us to compute the divergence (its
% transposed operator) and the Laplacian 
% (as the composition of the divergence and the gradient).

%%
% The gradient operator \(\nabla\) can be understood as a collection of 3 sparse 
% matrices \((\nabla_s)_{s=1,2,3}\) of size \((m,n)\)
% that computes each coordinate of \(\nabla u=(\nabla_s u)_{s=1,2,3}\)
% through the formula, for each face \(f\), 
% \[ (\nabla u)_f = \frac{1}{2A_f} \sum_{i \in f} u_i (N_f \wedge e_i) \]
% where \(A_f\) is the area of face \(f\), \(N_f\) is the normal to the
% face, \(e_i\) is the edge opposite to vertex \(i\), and
% \(\wedge\) is the cross product in \(\RR^3\).

%%
% Callback to get the coordinates of all the vertex of index \(i=1,2,3\) in
% all faces.

XF = @(i)X(:,F(i,:));

%%
% Compute un-normalized normal through the formula 
% \(e_1 \wedge e_2 \) where \(e_i\) are the edges.

Na = cross( XF(2)-XF(1), XF(3)-XF(1) );

%%
% Compute the area of each face as half the norm of the cross product.

amplitude = @(X)sqrt( sum( X.^2 ) );
A = amplitude(Na)/2;

%%
% Compute the set of unit-norm normals to each face.

normalize = @(X)X ./ repmat(amplitude(X), [3 1]);
N = normalize(Na);

%%
% Populate the sparse entries of the matrices for the 
% operator implementing \( \sum_{i \in f} u_i (N_f \wedge e_i) \).

I = []; J = []; V = []; % indexes to build the sparse matrices
for i=1:3
    % opposite edge e_i indexes
    s = mod(i,3)+1;
    t = mod(i+1,3)+1;
    % vector N_f^e_i
    wi = cross(XF(t)-XF(s),N);    
    % update the index listing
    I = [I, 1:m];
    J = [J, F(i,:)]; 
    V = [V, wi];
end

%% 
% Sparse matrix with entries \(1/(2A_f)\).

dA = spdiags(1./(2*A(:)),0,m,m); 

%%
% Compute gradient.

GradMat = {};
for k=1:3
    GradMat{k} = dA*sparse(I,J,V(k,:),m,n);
end

%%
% \(\nabla\) gradient operator.

Grad = @(u)[GradMat{1}*u, GradMat{2}*u, GradMat{3}*u]';

%%
% Compute divergence matrices as transposed of grad for the face area inner product.

dAf = spdiags(2*A(:),0,m,m); 
DivMat = {GradMat{1}'*dAf, GradMat{2}'*dAf, GradMat{3}'*dAf};

%%
% Div operator. 

Div = @(q)DivMat{1}*q(1,:)' + DivMat{2}*q(2,:)' + DivMat{3}*q(3,:)';

%%
% Laplacian operator as the composition of grad and div. 

Delta = DivMat{1}*GradMat{1} + DivMat{2}*GradMat{2} + DivMat{3}*GradMat{3};

%%
% Cotan of an angle between two vectors.

cota = @(a,b)cot( acos( dot(normalize(a),normalize(b)) ) );

%%
% Compute cotan weights Laplacian. 

I = []; J = []; V = []; % indexes to build the sparse matrices
Ia = []; Va = []; % area of vertices
for i=1:3
    % opposite edge e_i indexes
    s = mod(i,3)+1;
    t = mod(i+1,3)+1;
    % adjacent edge
    ctheta = cota(XF(s)-XF(i), XF(t)-XF(i));
    % ctheta = max(ctheta, 1e-2); % avoid degeneracy
    % update the index listing
    I = [I, F(s,:), F(t,:)];
    J = [J, F(t,:), F(s,:)]; 
    V = [V, ctheta, ctheta];
    % update the diagonal with area of face around vertices
    Ia = [Ia, F(i,:)];
    Va = [Va, A]; 
end
% Aread diagonal matrix
Ac = sparse(Ia,Ia,Va,n,n);
% Cotan weights
Wc = sparse(I,J,V,n,n);
% Laplacian with cotan weights.
DeltaCot = spdiags(full(sum(Wc))', 0, n,n) - Wc;

%%
% Check that the Laplacian with cotan weights is actually equal to 
% the composition of divergence and gradient.

fprintf('Should be 0: %e\n', norm(Delta-DeltaCot, 'fro')/norm(Delta, 'fro'));

%%
% Display a function \(f\) on the mesh.

f = X(2,:);
options.face_vertex_color = f(:);
clf; plot_mesh(X,F,options); 
axis('tight');
colormap parula(256);

%%
% Display its Laplacian.

g = Delta*f(:);
g = clamp(g, -3*std(g), 3*std(g));
options.face_vertex_color = rescale(g);
clf; plot_mesh(X,F,options); 
axis('tight');
colormap parula(256);

%% Heat Diffusion and Time Stepping
% The method developped in <#biblio [CraneWeischedelWardetzky13]> relies on
% the fact that the level set of the geodesic distance function to a
% starting point \(i\) agrees with the level set of the solution of
% the heat diffusion when the time of diffusion tends to zero. This
% fundamental result is proved in <#biblio [Varadhan67]>.

%%
% In fact, the
% same result holds true when replacing the heat diffusion solution by a
% single Euler implicit step in time, with time step \(t\). This means one
% should consider the solution \(u\) to the equation 
% \[ (\text{Id}+t \Delta) u = \delta_i \]
% where \(\delta_i\) is the Dirac vector at vertex index \(i\).

%%
% Select index \(i\).

i = 21000; 

%%
% Set time \(t\).

t = 1000;

%%
% Solve the linear system. 

delta = zeros(n,1); 
delta(i) = 1;
u = (Ac+t*Delta)\delta;

%%
% Display this solution.

options.face_vertex_color = u;
clf; plot_mesh(X,F,options); 
axis('tight');
colormap parula(256);

%EXO
%% Solve the heat diffusion equation 
%% \[ \frac{\partial g}{\partial t} = -\Delta g \]
%% by performing several implicit time stepping.
t = .1;
nlist = [1 3 10 50]; k = 1;
clf;
u = delta;
for i=1:max(nlist)
    % solve
    u = (Ac+t*DeltaCot)\u;
    % display
    if k<=4 && i==nlist(k)
        subplot(2,2,k); hold on;
        options.face_vertex_color = u;
        plot_mesh(X,F,options);
        axis('tight');
        colormap parula(256);
        title(['Iter ' num2str(i)]);
        drawnow;
        k = k+1;
    end
end
%EXO

%% Geodesic in Heat Method
% The main point of the method <#biblio [CraneWeischedelWardetzky13]> is to
% retrieve an approximation of the distance function \(\phi\) from the level set of
% implicit heat diffusion step \(u\). 

%%
% This is achieved by using the fact
% that \(\norm{\nabla \phi}=1\), i.e. one should have 
% \[ \nabla \phi \approx -\frac{\nabla u}{\norm{\nabla u}}. \]
% Solving this equation in the least square sense leads to the resolution
% of a Poisson equation. 

%%
% Compute the solution \(u\) with explicit time stepping.

t = .1;
u = (Ac+t*DeltaCot)\delta;

%%
% Compute the gradient field.

g = Grad(u); 

%%
% Normalize it to obtain 
% \[ h = -\frac{\nabla u}{\norm{\nabla u}}. \]

h = -normalize(g);

%%
% Integrate it back by solving
% \[ \Delta \phi = \text{div}(h).  \]

phi = Delta \ Div(h);

%%
% Display.

options.face_vertex_color = phi;
clf; plot_mesh(X,F,options); 
axis('tight');
colormap parula(256);

%% 
% Functions for displaying the level sets of the distance.

p = 30;
DispFunc = @(phi)cos(2*pi*p*phi);

%%
% Same, but using a different colormap.

options.face_vertex_color = DispFunc(phi);
clf; plot_mesh(X,F,options); 
axis('tight');
colormap parula(256);

%EXO
%% Display the result obtained for several time \(t\).
tlist = [.001 1 100 10000];
clf;
for k=1:length(tlist)
    t = tlist(k);
    % solve
    u = (Ac+t*DeltaCot)\delta;
    % re-integrate normalized gradient
    phi = -Delta \ Div( normalize( Grad(u) ) );
    % display
    options.face_vertex_color = DispFunc(phi);
    subplot(2,2,k);
    plot_mesh(X,F,options); 
    title(['t=' num2str(t)]);
    axis('tight');
    colormap parula(256);
end
%EXO

%EXO
%% Compute distances from an increasing number of starting points
%% that are computed using a farthest point sampling.
nlist = [1 2 3 4]; 
nlist = [3 6 10 20]; 
k = 1;
i = 21000; % for elephant (nose)
I = [i];
delta = zeros(n,1); 
delta(i) = 1;
t = 1/1000; 
clf;
for i=1:max(nlist)
    % solve
    u = (Ac+t*DeltaCot)\delta;
   	% re-integrate normalized gradient
    phi = -Delta \ Div( normalize( Grad(u) ) );
    % display
    if k<=4 && i==nlist(k)
        subplot(2,2,k); hold on;
        plot3(X(1,I), X(2,I), X(3,I), 'r.', 'MarkerSize', 25);
        options.face_vertex_color = DispFunc(phi);
        plot_mesh(X,F,options);
        drawnow;
        axis('tight');
        colormap parula(256);
        k = k+1;
    end
    % add fatherest point
    [~,j] = max(phi);
    delta(j) = 1;
    I(end+1) = j;
end
%EXO

%% Bibliography
% <html><a name="biblio"></a></html>

%%
% * [CraneWeischedelWardetzky13] K. Crane, C. Weischedel, M. Wardetzky, _Geodesics in heat: A new approach to computing distance based on heat flow_, ACM Transactions on Graphics , vol. 32, no. 5, pp. 152:1-152:11, 2013.
% * [Varadhan67] S.R.S. Varadhan, _On the behavior of the fundamental solution of the heat equation with variable coefficients_ Communications on Pure and Applied Mathematics 20, 2, 431-455, 1967

	



