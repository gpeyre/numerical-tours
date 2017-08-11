function plot_multiclasses(X,y,options)

% plot_multiclasses - display data for classification
%
%   plot_multiclasses(X,y,options);
%
%   Copyright (c) 2017 Gabriel Peyre

% helpers
ms = getoptions(options, 'ms', 8);
fs = getoptions(options, 'fs', 20);
disp_legend = getoptions(options, 'disp_legend', 1);
disp_dim = getoptions(options, 'disp_dim', 2);
edge_factor = getoptions(options, 'edge_factor', .25);
SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1], 'FontSize', fs);
Xm = @(X)X-repmat(mean(X,1), [size(X,1) 1]);
Cov = @(X)Xm(X)'*Xm(X);
col = [ [1 0 0]; [0 1 0]; [0 0 1]; [0 0 0]; [0 1 1]; [1 0 1]; [1 1 0]; ...
    [1 .5 .5]; [.5 1 .5]; [.5 .5 1]  ]';

[p,n] = size(X);
% list of classes
CL = unique(y);
k = length(CL);

% dimensionality reduction
if n>disp_dim
    [U,D,V] = svd(Xm(X),'econ');
    Z = Xm(X) * V;
else
    Z = X;
end

hold on;
lgd = {};
for i=1:min(k,size(col,2))
    I = find(y==CL(i));
    switch disp_dim
        case 2
            plot(Z(I,1), Z(I,2), 'o', 'MarkerFaceColor', col(:,i), 'MarkerEdgeColor', col(:,i)*edge_factor, 'MarkerSize', ms);
            % plot(Z(I,1), Z(I,2), '.', 'Color', col(:,i), 'MarkerSize', ms);
        case 3
            plot3(Z(I,1), Z(I,2), Z(I,3), 'o', 'MarkerFaceColor', col(:,i), 'MarkerEdgeColor', col(:,i)*edge_factor, 'MarkerSize', ms);
            % plot3(Z(I,1), Z(I,2), Z(I,3), '.', 'Color', col(:,i), 'MarkerSize', ms);
        otherwise
            error('Works only in 2D and 3D.');
    end
    lgd{end+1} = num2str(CL(i));
end
axis tight; axis equal; box on;
if disp_dim==3
    view(3); 
end
if disp_legend==1
    legend(lgd, 'Location', 'EastOutside');
end
SetAR(1);

end