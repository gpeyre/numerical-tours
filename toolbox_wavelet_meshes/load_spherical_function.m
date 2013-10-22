function f = load_spherical_function(name, pos, options)

% load_spherical_function - load a function on the sphere
%
%   f = load_spherical_function(name, pos, options);
%
%   Copyright (c) 2007 Gabriel Peyre

if iscell(pos)
    pos = pos{end};
end

if size(pos,1)>size(pos,2)
    pos = pos';
end

x = pos(1,:); x = x(:);

M = [];
if not(isstr(name))
    M = name; 
    name = 'image';
end

switch name
    case 'linear'        
        f = x;
    case 'cos'
        f = cos(6*pi*x);
    case 'singular'
        f = abs(x).^.4;
    case 'image'
        name = getoptions(options, 'image_name', 'lena');
        if isempty(M)
            M = load_image(name);
            q = size(M,1);
            q = min(q,512);
            M = crop(M,q);
            M = perform_blurring(M,4);
        end
        f = perform_spherical_interpolation(pos,M);
    case 'etopo'
        resol = 15;
        fname = ['ETOPO' num2str(resol)];
        fid = fopen(fname, 'rb');
        if fid<0
            error('Unable to read ETOPO file');
        end
        s = [360*(60/resol), 180*(60/resol)];     
        M = fread(fid, Inf, 'short');
        M = reshape(M, s(1),s(2) );
        fclose(fid);
        f = perform_spherical_interpolation(pos,M, 0);
    case {'earth' 'earth-grad'}
        filename = 'earth-bw';
        M = double( load_image(filename) );
        M = perform_blurring(M,4);
        if strcmp(name, 'earth-grad')
            G = grad(M);
            M = sqrt(sum(G.^2,3));
            M = perform_blurring(M,15);
        end
        f = perform_spherical_interpolation(pos,M,0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = perform_spherical_interpolation(pos,M,center)

if nargin<3
    center = 1;
end

qx = size(M,1);
qy = size(M,2);
Y = atan2(pos(2,:),pos(1,:))/(2*pi) + 1/2;
if center
    X = acos(pos(3,:))/(2*pi) + 1/4;
else
    X = acos(pos(3,:))/(pi);
end
x = linspace(0,1,qx);
y = linspace(0,1,qy);
f = interp2( y,x,M, Y(:),X(:) );