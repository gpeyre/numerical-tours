function [D,S,Q] = perform_fast_marching(W, start_points, options)

% perform_fast_marching - launch the Fast Marching algorithm, in 2D or 3D.
%
%   [D,S,Q] = perform_fast_marching(W, start_points, options)
%
%   W is an (n,n) (for 2D, d=2) or (n,n,n) (for 3D, d=3) 
%       weight matrix. The geodesics will follow regions where W is large.
%       W must be > 0.
%   'start_points' is a d x k array, start_points(:,i) is the ith starting point .
%
%   D is the distance function to the set of starting points.
%   S is the final state of the points : -1 for dead (ie the distance
%       has been computed), 0 for open (ie the distance is only a temporary
%       value), 1 for far (ie point not already computed). Distance function
%       for far points is Inf.
%   Q is the index of the closest point. Q is set to 0 for far points.
%       Q provide a Voronoi decomposition of the domain. 
%
%   Optional:
%   - You can provide special conditions for stop in options :
%       'options.end_points' : stop when these points are reached
%       'options.nb_iter_max' : stop when a given number of iterations is
%          reached.
%   - You can provide an heuristic in options.heuristic (typically that try to guess the distance
%       that remains from a given node to a given target).
%       This is an array of same size as W.
%   - You can provide a map L=options.constraint_map that reduce the set of
%       explored points. Only points with current distance smaller than L
%       will be expanded. Set some entries of L to -Inf to avoid any
%       exploration of these points.
%   - options.values set the initial distance value for starting points
%   (default value is 0).
%
%   See also: perform_fast_marching_3d.
%
%   Copyright (c) 2007 Gabriel Peyre


options.null = 0;

end_points = getoptions(options, 'end_points', []);
verbose = getoptions(options, 'verbose', 1);
nb_iter_max = getoptions(options, 'nb_iter_max', Inf);
values = getoptions(options, 'values', []);
L = getoptions(options, 'constraint_map', []);
H = getoptions(options, 'heuristic', []);
dmax = getoptions(options, 'dmax', Inf);

d = nb_dims(W);

if (d==4 && size(W,3)==2 && size(W,4)==2) || (d==4 && size(W,4)==6) || (d==5 && size(W,4)==3 && size(W,5)==3)
    % anisotropic fast marching
    if d==4 && size(W,3)==2 && size(W,4)==2
        % 2D vector field -> 3D field
        W1 = zeros(size(W,1), size(W,2), 3, 3);
        W1(:,:,1:2,1:2) = W; 
        W1(:,:,3,3) = 1;
        W = reshape(W1, [size(W,1) size(W,2), 1 3 3]);
        % convert to correct size
        W = cat(4, W(:,:,:,1,1), W(:,:,:,1,2), W(:,:,:,1,3), W(:,:,:,2,2), W(:,:,:,2,3), W(:,:,:,3,3) );        
    end
    if d==5
        % convert to correct size
        W = cat(4, W(:,:,:,1,1), W(:,:,:,1,2), W(:,:,:,1,3), W(:,:,:,2,2), W(:,:,:,2,3), W(:,:,:,3,3) );
    end
    
    if size(start_points,1)==2
        start_points(end+1,:) = 1;
    end
    if size(start_points,1)~=3
        error('start_points should be of size (3,n)');
    end
    
    % padd to avoid boundary problem
    W = cat(1, W(1,:,:,:), W, W(end,:,:,:));
    W = cat(2, W(:,1,:,:), W, W(:,end,:,:));
    W = cat(3, W(:,:,1,:), W, W(:,:,end,:));
    
%    if isempty(L)
        L = ones(size(W,1), size(W,2), size(W,3));
 %   end
 
    if dmax==Inf
        dmax = 1e15;
    end
    
%    start_points = start_points-1;
    alpha = 0;
    [D,Q] = perform_front_propagation_anisotropic(W, L, alpha, start_points,dmax);

    % remove boundary problems
    D = D(2:end-1,2:end-1,2:end-1);
    Q = Q(2:end-1,2:end-1,2:end-1);
    S = [];
    D(D>1e20) = Inf;
    return;
end


if d~=2 && d~=3
    error('Works only in 2D and 3D.');
end
if size(start_points,1)~=d
    error('start_points should be (d,k) dimensional with k=2 or 3.');
end

L(L==-Inf)=-1e9;
L(L==Inf)=1e9;
nb_iter_max = min(nb_iter_max, 1.2*max(size(W))^d);

% use fast C-coded version if possible
if d==2
    if exist('perform_front_propagation_2d')~=0
        [D,S,Q] = perform_front_propagation_2d(W,start_points-1,end_points-1,nb_iter_max, H, L, values);
    else
        error('You should compile the mex file, see compile_mex.m');
        [D,S] = perform_front_propagation_2d_slow(W,start_points,end_points,nb_iter_max, H);
        Q = [];
    end
elseif d==3
    [D,S,Q] = perform_front_propagation_3d(W,start_points-1,end_points-1,nb_iter_max, H, L, values);
end
Q = Q+1;

% replace C 'Inf' value (1e9) by Matlab Inf value.
D(D>1e8) = Inf;