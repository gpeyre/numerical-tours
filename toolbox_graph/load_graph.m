function [A,vertex] = load_graph(name,n,options)

% load_graph - load some graph
%
%   A is the adjacency sparse matrix.
%   vertex defines some nice embeding for gplot (watch out, it is of size 2 x nb_vertex).
%   name can be:
%       * '*.off' : load a mesh from a .off file
%       * 'square' : regular square
%       * 'cyclic' : cyclic graph
%       * 'random_network' a random network topology by the method
%           suggested by Waxman (1988):
%           - nodes are a Poisson process in the plane with scaled
%             Lebesgue mean measure 
%           - nodes u and v are connected with probability
%               P(u,v)=option.alpha*exp(-d(u, v)/(option.beta*L))
%             where alpha>0, beta<=1, d(u,v) is Euclidean distance, 
%             L is the maximum distance between any two nodes
%               option.alpha - maximal link probability
%               option.beta - parameter to control length of the edges. Increased <beta>
%                   yields a larger ratio of long edges to short edges
%
%   [A,vertex] = load_graph(name,n,options);
%
%   Copyright (c) 2005 Gabriel Peyré

options.null = 0;

if strcmp(name(end-2:end), 'off') || strcmp(name(end-2:end), 'ply') || strcmp(name(end-2:end), 'smf')
    % this is a mesh, load it
    [vertex,face] = read_mesh(name);
    vertex = vertex(1:2,:);
    A = sparse( triangulation2adjacency(face,vertex) );
    return;
end

switch lower(name)
    
    case 'perfect'
        
        t = linspace(0,2*pi, n+1); t(end) = [];
        vertex = [cos(t); sin(t)];
        A = ones(n) - eye(n);
    
    case 'square'
        if isfield(options, 'connectivity')
            connectivity = options.connectivity;
        else
            connectivity = 4;
        end
        [A,vertex] = gen_square_graph( round(sqrt(n)), connectivity );
        
    case 'cyclic'
        if isfield(options, 'generating')
            generating = options.generating;
        else
            generating = 4;
        end
        [A,vertex] = gen_cyclic_graph( n, generating );
        
        
    case 'chain'
        x = linspace(0,1,n);
        vertex = [x; x*0];
        i = [1:n 1:n]
        j = [1 1:n-1 2:n n];
        s =  i*0 + 1;
        A = sparse(i,j,s);
        
    case 'erdos'
        
        if isfield(options, 'proba')
            proba = options.proba;
        else
            proba = 0.1;
        end
        Kreg = 4;
        [A,vertex]=erdosRenyi(n,proba,Kreg);
        
    case 'rand'
        if isfield(options, 'proba')
            proba = options.proba;
        else
            proba = 0.1;
        end
        A = sparse( rand(n)<proba );
        A = max(A,A);
        % try to find some layout
		options.method = 'drawing'; options.A = A;
        vertex = perform_parameterization(0,0, options);
     
    case 'rand_triangulation'
        
        vertex = rand(2,n);
        face = delaunay(vertex(1,:), vertex(2,:));
        A = triangulation2adjacency(face,vertex);
    
    case 'randn_triangulation'
        
        vertex = randn(2,n);
        face = delaunay(vertex(1,:), vertex(2,:));
        A = triangulation2adjacency(face,vertex);
    
    case 'rand_clusters_triangulation'
        
        % centers
        s = 4.5;
        c = s * [0 0; 0.5 1; 1 0];
        vertex = [];
        for i=1:size(c,1)
            v = randn(2,round(n/3));
            v(1,:) = v(1,:) + c(i,1);
            v(2,:) = v(2,:) + c(i,2);
            vertex = [vertex v];
        end
        face = delaunay(vertex(1,:), vertex(2,:));
        A = triangulation2adjacency(face,vertex);
        
    case 'rand_network'

        if isfield(options, 'alpha')
            alpha = options.alpha;
        else
            alpha = 0.4;
        end
        if isfield(options, 'beta')
            beta = options.beta;
        else
            beta = 0.1;
        end
        [A, vertex]=waxtop(n, alpha, beta);
        
    otherwise
        
        [A,vertex] = load_graph_dataset(['data/' name]);
    
end

A = sparse(A);

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end




function [A,xy] = gen_cyclic_graph( n, S )

% gen_cyclic_graph - generate a Caley graph associated with
%   generating set S in Z/nZ.
%
%   [A,xy] = gen_cyclic_graph( n, S );
%
%   A : adjacency matrix of size nxn
%   xy : position (on the circle) of nx2.
%
%   Copyright (c) 2003 Gabriel Peyré

A = zeros(n,n);
x = (0:(n-1))'*2*pi/n;
xy = [ cos(x) sin(x) ];

for i=0:n-1
    for s = S
        A( i+1, mod(i+s,n)+1 ) = 1;
    end
end




function [A,xy] = gen_square_graph( n, connectivity )

% gen_square_graph - generate a Caley graph on a square grid, i.e.
%   of (Z/nZ)^2
%   associated to the generativ set S = { [-1,0], [1,0], [0,1], [0,-1] } 
%   for 'connectivity==4' and 
%   S = { [-1,-1], [1,1], [-1,1], [1,-1], [-1,0], [1,0], [0,1], [0,-1] } 
%   for 'connectivity==8'
%
%    [A,xy] = gen_square_graph( n, connectivity );
%
%   A : adjacency matrix of size n x n
%   xy : position matrix of size n x 2.
%
%   Copyright (c) 2003 Gabriel Peyré

if nargin<2
    connectivity = 4;
end

A = zeros(n^2,n^2);
xy = zeros(n^2,2);
h = 1/(n-1);

for i=0:n-1
for j=0:n-1
    k = i+n*j;
    xy(k+1,1) = i*h;
    xy(k+1,2) = j*h;
    if i<n-1
        A( k+1, k+2 ) = 1;
    end
    if i>0
        A( k+1, k ) = 1;
    end
    if j<n-1
        A( k+1, k+1+n ) = 1;
    end
    if j>0
        A( k+1, k+1-n ) = 1;
    end
    if connectivity==8
        if (i<n-1) & (j<n-1)
            A( k+1, k+2+n ) = 1;
        end
        if (i<n-1) & (j>0)
            A( k+1, k+2-n ) = 1;
        end
        if (i>0) & (j>0)
            A( k+1, k-n ) = 1;
        end
        if (i>0) & (j<n-1)
            A( k+1, k+n ) = 1;
        end
    end
end
end




function [A,xy] = load_graph_dataset(name)

% load_graph_dataset - load a graph data set
%
%   [A,xy] = load_graph_dataset(name);
%
%   name is the base name (eg. 'data/cal')
%   of two files 
%       [name '.cegde'] (data for edges)
%       [name '.cnode'] (data for nodes)
%
%   Copyright (c) 2005 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% node files

filename = [name, '.cnode'];
fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the node file.');
    return;
end

[xy,nbr_edges] = fscanf(fid,'%d %f %f');

fclose(fid);

xy = reshape(xy, 3, length(xy)/3);
xy = xy(2:3,:);
n = size(xy,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% edge file
filename = [name, '.cedge'];
fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the edge file.');
    return;
end

[C,nbr_edges] = fscanf(fid,'%d %d %d %f');

fclose(fid);

C = reshape(C, 4, length(C)/4);
d = C(end,:); % L^2 dist
C = C(2:3,:) + 1;

% number of edges
if max(C(:))>n
    error('Edge contains index larger than number of vertices');
end
p = size(C,2);

A = spalloc(n,n,p);
for i=1:p
    A(C(1,i),C(2,i)) = d(i);
end


function [A,vertex]=erdosRenyi(nv,p,Kreg)
%Funciton [G]=edosRenyi(nv,p,Kreg) generates a random graph based on
%the Erdos and Renyi algoritm where all possible pairs of 'nv' nodes are
%connected with probability 'p'. 
%
% Inputs:
%   nv - number of nodes 
%   p  - rewiring probability
%   Kreg - initial node degree of for regular graph (use 1 or even numbers)
%
% Output:
%   G is a structure inplemented as data structure in this as well as other
%   graph theory algorithms.
%   G.Adj   - is the adjacency matrix (1 for connected nodes, 0 otherwise).
%   G.x and G.y -   are row vectors of size nv wiht the (x,y) coordinates of
%                   each node of G.
%   G.nv    - number of vertices in G
%   G.ne    - number of edges in G
%
%Created by Pablo Blinder. blinderp@bgu.ac.il
%
%Last update 25/01/2005

%build regular lattice 
A=sparse(nv,nv);
Kreg=fix(abs(Kreg)/2);Kreg=(Kreg<1)+Kreg;

for k=1:Kreg
    A=sparse(A+diag(ones(1,length(diag(A,k))),k)+diag(ones(1,length(diag(A,nv-k))),nv-k));
end
ne0=nnz(A);
%find connected pairs
[v1,v2]=find(A);
% P=permPairs(nv);%my version is faster
Dis=(rand(length(v1),1)<=p);%pairs to disconnect
A(v1(Dis),v2(Dis))=0;
vDis=unique([v1(Dis),v2(Dis)]);%disconnected vertices
nDis=ne0-nnz(A);sum(Dis);

%cycle trough disconnected pairs
disconPairs=[v1(Dis),v2(Dis)];
for n=1:nDis
    %choose one of the vertices from the disconnected pair
    i=ceil(rand*size(disconPairs,1));
    j=logical(1+rand>0.5);
    vDisToRec=disconPairs(i,j);
    %find non adjacent vertices and reconnect
    adj=[find(A(:,vDisToRec)) ; find(A(vDisToRec,:))'];
    nonAdj=setdiff(1:nv,adj);
    vToRec=nonAdj(ceil(rand*length(nonAdj)));
    S=sort([vDisToRec vToRec]);
    A(S(1),S(2))=1;
end
[x,y]=getNodeCoordinates(nv);
%make adjacency matrix symetric
A=A+fliplr((flipud(triu(A))));
% G=struct('Adj',A,'x',x','y',y','nv',nv,'ne',nnz(A));
vertex = [x,y];

function [x,y]=getNodeCoordinates(nv)
%Adapted from circle.m by Zhenhai Wang <zhenhai@ieee.org>. For more details
%see under  MATLAB Central >  File Exchange > Graphics > Specialized
%Plot and Graph Types > Draw a circle.

center=[0,0];
theta=linspace(0,2*pi,nv+1);
rho=ones(1,nv+1);%fit radius and nv
[X,Y] = pol2cart(theta',rho');
X=X+center(1);
Y=Y+center(2);
x=X(1:end-1)*10;
y=Y(1:end-1)*10;


function P=permPairs(N)
%Produces all pairs of pairs from 1 to N.
%It is ~30% to 50% faster that nchoosek(1:N,2).
%Created by Pablo 02/12/2003

ini_i=1;
ini_j=ini_i+1;
r0=1;
P=[];
for i=ini_i:N-1
    lj=N-ini_i;
    P(:,r0:lj+r0-1)=[ones(1,lj)*i;ini_j:N];
    r0=r0+lj;
    ini_i=ini_i+1;
    ini_j=ini_i+1;
end
P=P';






function [adj_matr, n_coord]=waxtop(n_points, alpha, beta, domain) 
% WAXTOP Simulate and plot a random network topology by the method
%   suggested by Waxman (1988):
%   - nodes are a Poisson process in the plane with scaled
%    Lebesgue mean measure 
%   - nodes u and v are connected with probability
%    P(u,v)=alpha*exp(-d(u, v)/(beta*L)),
%    where alpha>0, beta<=1, d(u,v) is Euclidean distance, 
%    L is the maximum distance between any two nodes
%
% [A, vertex]=waxtop(lambda, alpha, beta, domain);
%
% Inputs: 
%   lambda - intensity of the Poisson process
%   alpha - maximal link probability
%   beta - parameter to control length of the edges. Increased <beta>
%     yields a larger ratio of long edges to short edges
%   domain - bounds for the region. A 4-dimensional vector in
%     the form [x_min x_max y_min y_max]. 
%
% Outputs:
%   adj_matr - adjacency matrix of the graph of the topology
%   n_coord - coordinates of the nodes

% Authors: R.Gaigalas, I.Kaj
% v1.4 07-Nov-01

if nargin<1 % default parameter values
    lambda = 0.6; % intensity of the Poisson process
end
if nargin<2
    alpha = 0.4; % parameter for the link probability
end
if nargin<3
    beta = 0.1; % parameter for the link probability
end
if nargin<4
    domain = [0 10 0 10]; % bounds for the "geografical" domain
end

x_min = domain(1);
x_max = domain(2);
y_min = domain(3);
y_max = domain(4);
clear domain;

% number of points is Poisson distributed
% with intensity proportional to the area
area = (x_max-x_min)*(y_max-y_min);
% n_points = poissrnd(lambda*area)


% given the number of points, nodes are uniformly distributed
n_coord = rand(n_points, 2);
n_coord(:,1) = n_coord(:,1)*(x_max-x_min)+x_min;
n_coord(:,2) = n_coord(:,2)*(y_max-y_min)+y_min;

% create a matrix with distances
% x_rep = n_coord(:, 1)*ones(1, n_points);
% y_rep = n_coord(:, 2)*ones(1, n_points);
x_rep = repmat(n_coord(:, 1), 1, n_points);
y_rep = repmat(n_coord(:, 2), 1, n_points);
dist_matr = sparse(triu(((x_rep-x_rep').^2+(y_rep-y_rep').^2).^0.5, 1));

% create the matrix of probabilities
prob_matr = alpha*spfun('exp', -dist_matr./(beta*max(max(dist_matr))));

% generate the adjacency matrix
coin = sprand(dist_matr);
adj_matr = (coin>0) & (coin < prob_matr);

% test for connectivity
% s_matr = speye(size(adj_matr));
% for i=1:n_points-1
%   s_matr = s_matr+adj_matr^i;
% end
% length(find(s_matr==0))

% plot the network
plot(n_coord(:,1), n_coord(:,2), '.');

hold on;
gplot(adj_matr, n_coord);
hold off;
