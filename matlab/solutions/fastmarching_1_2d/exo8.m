x1 = [ [175;5] [21;42] [48;133] [244;78] [191;40] ...
         [100;13] [66;42] [183;66] [220;117]];
% gradient
D1 = D; D1(D1==Inf) = max(D1(D1~=Inf));
G = grad(D1, options);
G = G ./ repmat( sqrt( sum(G.^2, 3)+1e-9 ), [1 1 2]);
% extract centerlines
curves = {};
for k=1:size(x1,2)
    % extract curve
    gamma = x1(:,k);
    for iter=1:1.5*n/tau
        g = [interp2(1:n,1:n,G(:,:,1),gamma(2,end),gamma(1,end)); ...
            interp2(1:n,1:n,G(:,:,2),gamma(2,end),gamma(1,end)) ];
        gamma(:,end+1) = clamp( gamma(:,end) - tau*g, 1,n );
        if norm(gamma(:,end)-x0)<1
            break;
        end
    end
    gamma(:,end+1) = x0;
    curves{end+1} = gamma;
end
% display the curves
clf;
hold on;
imageplot(f, 'Boundary points');
for i=1:length(curves)
    c = curves{i};
    h = plot(c(2,:),c(1,:), 'r'); set(h, 'LineWidth', 2);    
end
h = plot(x0(2,:),x0(1,:), '.g'); set(h, 'MarkerSize', 25);
h = plot(x1(2,:), x1(1,:), '.b'); set(h, 'MarkerSize', 25);
