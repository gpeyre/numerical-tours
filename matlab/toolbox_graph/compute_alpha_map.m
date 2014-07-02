function [a,c] = compute_alpha_map(aname, options)

% compute_alpha_map - set up alphamap
%
%   [a,c] = compute_alpha_map(aname, options);
%
%   aname is the name of the alphamap can be 'rampup', 'rampdown', 'gaussian'.
%
%   a is an alphamap.
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<2
    cname = 'jet';
end

options.null = 0;
nalpha = getoptions(options, 'nalpha',256);
sigma = getoptions(options, 'sigma',.08);
c = getoptions(options, 'center',.4);

x = linspace(0,1,nalpha);

switch lower(aname)
    case 'rampup'
        a = x;
    case 'rampdown'
        a = 1-x;
    case 'gaussian'
        a = exp( -(x-c).^2 / (2*sigma^2) );
end

alphamap(a);

cname = getoptions(options, 'cname', 'jet');
% compute an associate colormap 
c0 = eval([cname '(256);']);
nc = size(c0,1);
si = rescale(cumsum(a));
s = linspace(0,1,nc);
c = zeros(length(si),3);
for i=1:3
    c(:,i) = interp1(s,c0(:,i),si);
end