% \(phi\).
gD = grad(phi,options);
d = max(eps, sqrt(sum(gD.^2,3)) );
g = gD ./ repmat( d, [1 1 2] );
% gradient
G = d .* div( g,options ) - lambda*(f0-c1).^2 + lambda*(f0-c2).^2;
