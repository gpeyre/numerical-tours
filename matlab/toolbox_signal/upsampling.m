function y = upsampling(x,d,p)

% upsampling - add p zeros between samples along dimension d
%
%   y = upsampling(x,d,p);
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
        y = zeros(p*size(x,1),size(x,2),size(x,3)); 
        y(1:p:end,:,:) = x;
    case 2
        y = zeros(size(x,1),p*size(x,2),size(x,3)); 
        y(:,1:p:end,:) = x;
    case 3
        y = zeros(size(x,1),size(x,2),p*size(x,3)); 
        y(:,:,1:p:end) = x;
    otherwise
        error('Not implemented');
end
