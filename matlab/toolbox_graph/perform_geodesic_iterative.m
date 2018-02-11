function [U,err,Usvg] = perform_geodesic_iterative(vertex, faces, W, I, options)

% perform_geodesic_iterative - compute the geodesic on a mesh using an iterative scheme
%
%   [U,err,Usvg] = perform_geodesic_iterative(vertex, faces, W, I, options);
%
%   INPUTS:
%   vertex and faces describes a mesh in arbitrary dimension.
%   W(i) is the weight (metric) at vertex indexed by i.
%   I is a list of index of starting points.
%   options.niter is the number of iteration (early quit if 
%       convergence is reached)
%   options.verb=1 to show progression
%   options.U gives an initialization. 
%   options.svg_rate gives the rate at wich Usvg is filled.
%
%   OUTPUT:
%   U(i) is the geodesic distance between point of index i and I.
%   err(k) is norm(U_{k+1}-U_k) where U_k is the solution at iteration k.
%       err should converge to zero.
%   Usvg(:,k) is the solution at iteration options.svg_rate*k.
%
%   Warning: to ensure convergence, options.U should be smaller than the
%   final solution (monotone convergence).
%
%   Copyright (c) 2011 Gabriel Peyre


options.null = 0;
verb = getoptions(options, 'verb', 1);
svg_rate = getoptions(options, 'svg_rate', 10);
niter = getoptions(options, 'niter', 200);

dotp = @(u,v)sum(u.*v,1);
R = @(u)reshape(u, [1 1 length(u)]);
Inv1 = @(M,d)[M(2,2,:)./d -M(1,2,:)./d; -M(2,1,:)./d M(1,1,:)./d];
Inv  = @(M)Inv1(M, M(1,1,:).*M(2,2,:) - M(1,2,:).*M(2,1,:));
Mult = @(M,u)[M(1,1,:).*u(1,1,:) + M(1,2,:).*u(2,1,:);  M(2,1,:).*u(1,1,:) + M(2,2,:).*u(2,1,:)];

n = size(vertex,2);

i = [faces(1,:) faces(2,:) faces(3,:) ];
j = [faces(2,:) faces(3,:) faces(1,:) ];
k = [faces(3,:) faces(1,:) faces(2,:) ];

err = [];
U = getoptions(options, 'U', []);
if isempty(U)
    U = zeros(n,1);
end
Usvg = [];


x  = vertex(:,i);
x1 = vertex(:,j) - x;
x2 = vertex(:,k) - x;
% inner product matrix
C = [R(dotp(x1,x1)) R(dotp(x1,x2)); ...
    R(dotp(x2,x1)) R(dotp(x2,x2))];
S = Inv(C);

% a = <S*1,1>
a = sum(sum(S));

w = R( W(i) );

% edge length
L1 = sqrt(dotp(x1,x1)); L1 = L1(:).*w(:); 
L2 = sqrt(dotp(x2,x2)); L2 = L2(:).*w(:); 


for it=1:niter
    if verb
        progressbar(it,niter);
    end
    uj = U(j);
    uk = U(k);
    u = [R(uj); R(uk)];
    
    % b = <S*1,u>
    b = dotp( sum(S,2), u );
    % c = <S*u,u> - W.^2;
    c = dotp( Mult(S,u), u ) - w.^2;
    % delta = b^2 - a*c
    delta = max( b.^2 - a.*c, 0);
    % solution
    d = (b + sqrt(delta) )./a;    
    
    % direction of the update
    % g=X*alpha,  alpha = S*(u-d*1)
    alpha = Mult( S, u - repmat(d, 2, 1) );
    J = find( alpha(1,1,:)>0 | alpha(2,1,:)>0 );    
    
    % update along edges
    
    d1 = L1 + uj(:);
    d2 = L2 + uk(:);
    d = d(:); 
    d(J) = min(d1(J), d2(J));
    
    U1 = accumarray(i', d, [n 1], @min);  U1(U1==0) = Inf;
    % boundary condition
    U1(I) = 0;
    % enforce monotony
    if min(U1-U)<-1e-5
      %  warning('Monotony problem');
    end
    %
    err(end+1) = norm(U-U1, 'fro');
    if err(end)==0
        break;
    end
    % update
    U = U1;
    if mod(it,svg_rate)==1 && nargout>2
        Usvg(:,end+1) = U;
    end
end