function [D,Dsvg,Ssvg] = perform_dijstra_fm(W, pstart, options)

% perform_fm_dijstra - *slow* (matlab) implementation of Dijstra and FM
%
%   [D,Dsvg,Ssvg] = perform_fm_dijstra(W, pstart, options);
%
%   W is an (n,n) metric matrix.
%   pstart is a (2,k) starting points.
%   options.method is either 'fm' or 'dijstra'
%   
%   D is the final distance map to pstart
%   options.svg_rate gives the rate at wich Dsvg and Ssvg is filled.
%   options.niter can be used to limit the total number of steps (partial propagation). 
%   
%   Copyright (c) 2012 Gabriel Peyre

method = getoptions(options, 'method', 'fm');
svg_rate = getoptions(options, 'svg_rate', 10);
niter = getoptions(options, 'niter', Inf);
bound = getoptions(options, 'bound', 'sym');

%%
% Size.

n = size(W,1);

%%
% The four displacement vector to go to the four neightbors.

neigh = [[1;0] [-1;0] [0;1] [0;-1]];

%%
% For simplicity of implementation, we use periodic boundary conditions.

switch bound
    case 'per'
        boundary = @(x)mod(x-1,n)+1;
    case 'sym'
        boundary = @(x)x.*(x<=n & x>0) + (2-x).*(x<=0) + (2*n-x).*(x>n);
    otherwise
        error('Works only for per and sym.');
end

%%
% For a given grid index |k|, and a given neighboring index k in \({1,2,3,4}\), 
% |Neigh(k,i)| gives the corresponding grid neigboring index.

ind2sub1 = @(k)[rem(k-1, n)+1; (k - rem(k-1, n) - 1)/n + 1]; 
sub2ind1 = @(u)(u(2)-1)*n+u(1);
Neigh = @(k,i)sub2ind1( boundary(ind2sub1(k)+neigh(:,i)) );

%%
% Stack of starting points.

for i=1:size(pstart,2)
    I(i) = sub2ind1(pstart(:,i));
end


%%
% Initialize the distance to \(+\infty\), excepted for the boundary conditions.

D = zeros(n)+Inf; % current distance
D(I) = 0; 


%%
% Initialize the state to 0 (unexplored), excepted for the boundary point to \(1\)
% (front).

S = zeros(n);
S(I) = 1; % open

%%
% Run!

iter = 0;
Dsvg = []; Ssvg = [];
while not(isempty(I)) && iter<=niter
    iter = iter+1;
    % pop from stack
    [tmp,j] = sort(D(I)); j = j(1);
    i = I(j); I(j) = [];
    % declare dead
    S(i) = -1; 
    % The list of the four neighbors
    J = [Neigh(i,1); Neigh(i,2); Neigh(i,3); Neigh(i,4)];
    % Remove those that are dead
    J(S(J)==-1) = [];
    % Add them to the front
    J1 = J(S(J)==0);
    I = [I; J1];
    S(J1) = 1;
    % update neighbor values
    for j=J'
        dx = min( D([Neigh(j,1) Neigh(j,2)]) );
        dy = min( D([Neigh(j,3) Neigh(j,4)]) ); 
        switch method
            case 'dijkstra'
                % Dijkstra update
                D(j) = min(dx+W(j), dy+W(j));
            case 'fm'
                % FM update                
                Delta = 2*W(j) - (dx-dy)^2;
                if Delta>=0
                    D(j) = (dx+dy+sqrt(Delta))/2;
                else
                    D(j) = min(dx+W(j), dy+W(j));
                end
            otherwise
                error('Only dijstra and fm are allowed.');
        end
    end
    % svd
    if mod(iter,svg_rate)==1 && nargout>2
        Dsvg(:,:,end+1) = D;
        Ssvg(:,:,end+1) = S;
    end
end