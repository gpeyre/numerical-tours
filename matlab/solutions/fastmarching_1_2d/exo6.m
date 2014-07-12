options.nb_iter_max = Inf;
options.end_points = [];
tau = .8;
curves = {};
for i=1:2
    % FM
    [D,S,Q] = perform_fast_marching(1./W, x0(:,i), options);
    % gradient
%    D1 = D; D1(D1==Inf) = max(D1(D1~=Inf));
    G = grad(D, options);
    G = G ./ repmat( sqrt( sum(G.^2, 3)+1e-9 ), [1 1 2]);
    % Geodesic
    for j=1:2
        gamma = x1(:,j);
        for iter=1:1.5*n/tau
            g = [interp2(1:n,1:n,G(:,:,1),gamma(2,end),gamma(1,end)); ...
                interp2(1:n,1:n,G(:,:,2),gamma(2,end),gamma(1,end)) ];
            gamma(:,end+1) = gamma(:,end) - tau*g;
            if norm(gamma(:,end)-x0(:,i))<1
                break;
            end
        end
        gamma(:,end+1) = x0(:,i);
        curves{end+1} = gamma;
    end    
end
% display the curves
col = {'r', 'g', 'b', 'c'};
clf;
hold on;
imageplot(A, 'Boundary points');
for i=1:length(curves)
    c = curves{i};
    h = plot(c(2,:),c(1,:), col{i}); set(h, 'LineWidth', 2);    
end
h = plot(x0(2,:),x0(1,:), '.g'); set(h, 'MarkerSize', 25);
h = plot(x1(2,:), x1(1,:), '.b'); set(h, 'MarkerSize', 25);
