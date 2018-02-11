gD = grad(phi,options);
% normalized gradient
d = max(eps, sqrt(sum(gD.^2,3)) );
g = gD ./ repmat( d, [1 1 2] );
% gradient
G = - W .* d .* div( g,options ) - sum(gW.*gD,3);
