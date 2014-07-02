function y = subsampling(x,d,p)

% downsampling - subsampling along dimension d
%
%   y = subsampling(x,d,p);
%
%   default is p==2, d==1
%
%   Copyright (c) 2009 Gabriel Peyre

if nargin<3
    p = 2;
end
if nargin<2
    d = 1;
end

switch d
    case 1
        y = x(1:p:end,:,:);
    case 2
        y = x(:,1:p:end,:);
    case 3
        y = x(:,:,1:p:end);
    otherwise
        error('Not implemented');
end
