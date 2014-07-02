function TV = compute_total_variation(y, options)

% compute_total_variation - compute the total variation of an image
%
%   TV = compute_total_variation(y, options);
%
%   See also: perform_tv_projection, grad, div.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;

nbdims = 2;
if size(y,1)==1 || size(y,2)==1
    nbdims = 1;
end
if size(y,1)>1 && size(y,2)>1 && size(y,3)>1
    nbdims = 3;
end

if nbdims==1
    TV = sum( abs(diff(y)) );
    return;
end

weight_tv = getoptions(options, 'weight_tv', ones(size(y)) );
tv_norm = getoptions(options,'tv_norm', 'l2');

% options.bound = 'per';
g = grad(y, options);

switch tv_norm
    case 'l2'
        TV = sqrt( sum(g.^2,nbdims+1) );
    case 'l1'
        TV = sum( abs(g),nbdims+1);
    case 'linf'
        TV = max( abs(g),[],nbdims+1);
    otherwise
        error('Unknown norm.');
end
TV = TV .* weight_tv;
TV = sum(TV(:));