function [M,Normal] = load_gim(name, options)

% load_gim - load a geometry image, either from a file or synthetic.
%
%   [M,Normal] = load_gim(name, rep, type);
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if isfield(options, 'rep')
    rep = options.rep;
else    
    rep = '';
end
if isfield(options, 'type')
    type = options.type;
else
    type = 'gim';
end
if isfield(options, 'n')
    n = options.n;
else
    n = 256;
end

switch lower(name)

    case 'torus'

        if isfield(options, 'a')
            a = options.a;
        else
            a = 0.3;
        end
        c = 1-a;
        x = linspace(0,2*pi,n);
        [Y,X] = meshgrid(x,x);
        M = zeros(n,n,3);
        M(:,:,1) = 1+(c+a*cos(X)).*cos(Y);
        M(:,:,2) = 1+(c+a*cos(X)).*sin(Y);
        M(:,:,3) = a + a*sin(X);
        % tangent
        Tx = zeros(n,n,3);
        Tx(:,:,1) = -a*sin(X).*cos(Y);
        Tx(:,:,2) = -a*sin(X).*sin(Y);
        Tx(:,:,3) = a*cos(X);
        Ty = zeros(n,n,3);
        Ty(:,:,1) = -(c+a*cos(X)).*sin(Y);
        Ty(:,:,2) = (c+a*cos(X)).*cos(Y);
        Ty(:,:,3) = 0;
        % normal
        Normal = zeros(n,n,3);
        Normal(:,:,1) = Tx(:,:,2).*Ty(:,:,3)-Tx(:,:,3).*Ty(:,:,2);
        Normal(:,:,2) = -Tx(:,:,1).*Ty(:,:,3)+Tx(:,:,3).*Ty(:,:,1);
        Normal(:,:,3) = Tx(:,:,1).*Ty(:,:,2)-Tx(:,:,2).*Ty(:,:,1);
        Normal = Normal ./ repmat( sqrt( sum(Normal.^2,3) ), [1 1 3]);

    otherwise

        filename = [rep name '-' type '.gim'];
        M = read_gim(filename);
        filename = [rep name '-' type '.normal.gim'];
        Normal = read_gim(filename);
        Normal = Normal*2-1;
        n = size(M,1);

end