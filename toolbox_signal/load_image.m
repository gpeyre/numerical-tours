function M = load_image(type, n, options)

% load_image - load benchmark images.
%
%   M = load_image(name, n, options);
%
%   name can be:
%   Synthetic images:
%       'chessboard1', 'chessboard', 'square', 'squareregular', 'disk', 'diskregular', 'quaterdisk', '3contours', 'line',
%       'line_vertical', 'line_horizontal', 'line_diagonal', 'line_circle',
%       'parabola', 'sin', 'phantom', 'circ_oscil',
%       'fnoise' (1/f^alpha noise).
%   Natural images:
%       'boat', 'lena', 'goldhill', 'mandrill', 'maurice', 'polygons_blurred', or your own.
%   
%   Copyright (c) 2004 Gabriel Peyre

if nargin<2
    n = 512;
end
options.null = 0;

if iscell(type)
    for i=1:length(type)
        M{i} = load_image(type{i},n,options);
    end
    return;
end

type = lower(type);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for geometric objects
eta         = getoptions(options, 'eta', .1);
gamma       = getoptions(options, 'gamma', 1/sqrt(2));
radius      = getoptions(options, 'radius', 10);
center      = getoptions(options, 'center', [0 0]);
center1     = getoptions(options, 'center1', [0 0]);
w           = getoptions(options, 'tube_width', 0.06);
nb_points   = getoptions(options, 'nb_points', 9);
scaling     = getoptions(options, 'scaling', 1);
theta       = getoptions(options, 'theta', 30 * 2*pi/360);
eccentricity = getoptions(options, 'eccentricity', 1.3);
sigma = getoptions(options, 'sigma', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for the line, can be vertical / horizontal / diagonal / any
if strcmp(type, 'line_vertical')
    eta = 0.5;              % translation
    gamma = 0;      % slope
elseif strcmp(type, 'line_horizontal')
    eta = 0.5;              % translation
    gamma = Inf;      % slope
elseif strcmp(type, 'line_diagonal')
    eta = 0;              % translation
    gamma = 1;      % slope
end

if strcmp(type(1:min(12,end)), 'square-tube-')
    k = str2double(type(13:end));
    c1 = [.22 .5]; c2 = [1-c1(1) .5];
    eta = 1.5;
    r1 = [c1 c1] + .21*[-1 -eta 1 eta];
    r2 = [c2 c2] + .21*[-1 -eta 1 eta];
    M = double( draw_rectangle(r1,n) | draw_rectangle(r2,n) );
    if mod(k,2)==0
        sel = n/2-k/2+1:n/2+k/2;
    else
        sel = n/2-(k-1)/2:n/2+(k-1)/2;        
    end
    M( round(.25*n:.75*n), sel ) = 1;
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch lower(type)
    
    case 'constant'
        M = ones(n);
    
    case 'ramp'
        x = linspace(0,1,n);
        [Y,M] = meshgrid(x,x);
        
    case 'bump'
        
        s = getoptions(options, 'bump_size', .5);
        c = getoptions(options, 'center', [0 0]);
        if length(s)==1
            s = [s s];
        end
        x = linspace(-1,1,n);
        [Y,X] = meshgrid(x,x);
        X = (X-c(1))/s(1); Y = (Y-c(2))/s(2);
        M = exp( -(X.^2+Y.^2)/2 );
        
    case 'periodic'
        x = linspace(-pi,pi,n)/1.1;
        [Y,X] = meshgrid(x,x);
        f = getoptions(options, 'freq', 6);
        M = (1+cos(f*X)).*(1+cos(f*Y));
        
    case {'letter-x' 'letter-v' 'letter-z' 'letter-y'}
        M = create_letter(type(8), radius, n);
        
    case 'l'
        r1 = [.1 .1  .3 .9];
        r2 = [.1 .1 .9 .4];
        M = double( draw_rectangle(r1,n) | draw_rectangle(r2,n) );
        
    case 'ellipse'
        c1 = [0.15 0.5];
        c2 = [0.85 0.5];
        x = linspace(0,1,n);
        [Y,X] = meshgrid(x,x);
        d = sqrt((X-c1(1)).^2 + (Y-c1(2)).^2) + sqrt((X-c2(1)).^2 + (Y-c2(2)).^2);
        M = double( d<=eccentricity*sqrt( sum((c1-c2).^2) ) );
    case 'ellipse-thin'
        options.eccentricity = 1.06;
        M = load_image('ellipse', n, options);
    case 'ellipse-fat'
        options.eccentricity = 1.3;
        M = load_image('ellipse', n, options);
        
    case 'square-tube'
        c1 = [.25 .5];
        c2 = [.75 .5];
        r1 = [c1 c1] + .18*[-1 -1 1 1];
        r2 = [c2 c2] + .18*[-1 -1 1 1];        
        r3 = [c1(1)-w c1(2)-w c2(1)+w c2(2)+w];
        M = double( draw_rectangle(r1,n) | draw_rectangle(r2,n) | draw_rectangle(r3,n) );
    case 'square-tube-1'
        options.tube_width = 0.03;
        M = load_image('square-tube', n, options);        
    case 'square-tube-2'
        options.tube_width = 0.06;
        M = load_image('square-tube', n, options);
    case 'square-tube-3'
        options.im = 0.09;
        M = load_image('square-tube', n, options);
    case 'polygon'
        theta = sort( rand(nb_points,1)*2*pi );
        radius = scaling*rescale(rand(nb_points,1), 0.1, 0.93);        
        points = [cos(theta) sin(theta)] .* repmat(radius, 1,2);
        points = (points+1)/2*(n-1)+1; points(end+1,:) = points(1,:);        
        M = draw_polygons(zeros(n),0.8,{points'});
        [x,y] = ind2sub(size(M),find(M));
        p = 100; m = length(x);
        lambda = linspace(0,1,p);
        X = n/2 + repmat(x-n/2, [1 p]) .* repmat(lambda, [m 1]);
        Y = n/2 + repmat(y-n/2, [1 p]) .* repmat(lambda, [m 1]);
        I = round(X) + (round(Y)-1)*n;
        M = zeros(n); M(I) = 1;
    case 'polygon-8'
        options.nb_points = 8;
        M = load_image('polygon', n, options);
    case 'polygon-10'
        options.nb_points = 8;
        M = load_image('polygon', n, options);
    case 'polygon-12'
        options.nb_points = 8;
        M = load_image('polygon', n, options);
    case 'pacman'
        options.radius = 0.45;
        options.center = [.5 .5];
        M = load_image('disk', n, options);        
        x = linspace(-1,1,n);
        [Y,X] = meshgrid(x,x);
        T =atan2(Y,X);
        M = M .* (1-(abs(T-pi/2)<theta/2));
    case 'square-hole'
        options.radius = 0.45;
        M = load_image('disk', n, options); 
        options.scaling = 0.5;
        M = M - load_image('polygon-10', n, options); 
        
    case 'grid-circles'
        if isempty(n)
            n = 256;
        end
        f = getoptions(options, 'frequency', 30);
        eta = getoptions(options, 'width', .3);
        x = linspace(-n/2,n/2,n) - round(n*0.03);
        y = linspace(0,n,n);
        [Y,X] = meshgrid(y,x);
        R = sqrt(X.^2+Y.^2);
        theta = 0.05*pi/2;
        X1 = cos(theta)*X+sin(theta)*Y;
        Y1 = -sin(theta)*X+cos(theta)*Y;
        A1 = abs(cos(2*pi*R/f))<eta;
        A2 = max( abs(cos(2*pi*X1/f))<eta, abs(cos(2*pi*Y1/f))<eta );

        M = A1;
        M(X1>0) = A2(X1>0);
        
    case 'chessboard1'
        x = -1:2/(n-1):1;
        [Y,X] = meshgrid(x,x);
        M = (2*(Y>=0)-1).*(2*(X>=0)-1);
        
    case 'chessboard'
        width = getoptions(options, 'width', round(n/16) );
        [Y,X] = meshgrid(0:n-1,0:n-1);
        M = mod( floor(X/width)+floor(Y/width), 2 ) == 0;
        
    case 'square'
        if ~isfield( options, 'radius' )
            radius = 0.6;
        end
        x = linspace(-1,1,n);
        [Y,X] = meshgrid(x,x);
        M = max( abs(X),abs(Y) )<radius;
        
    case 'squareregular'
        M = rescale(load_image('square',n,options));
        if not(isfield(options, 'alpha'))
            options.alpha = 3;
        end
        S = load_image('fnoise',n,options);
        M = M + rescale(S,-0.3,0.3); 
        
    case 'regular1'
        options.alpha = 1;
        M = load_image('fnoise',n,options);
    case 'regular2'
        options.alpha = 2;
        M = load_image('fnoise',n,options);
    case 'regular3'
        options.alpha = 3;
        M = load_image('fnoise',n,options);
        
    case 'sparsecurves'
        options.alpha = 3;
        M = load_image('fnoise',n,options);
        M = rescale(M);
        ncurves = 3;
        M = cos(2*pi*ncurves);
        
    case 'geometrical'

        J = getoptions(options, 'Jgeometrical', 4);
        sgeom = 100*n/256;
        options.bound = 'per';
        A = ones(n);
        for j=0:J-1
            B = A;
            for k=1:2^j
                I = find(B==k);
                U = perform_blurring(randn(n),sgeom,options);
                s = median(U(I));
                I1 = find( (B==k) & (U>s) );
                I2 = find( (B==k) & (U<=s) );
                A(I1) = 2*k-1;
                A(I2) = 2*k;
            end
        end
        M = A;
        
    case 'lic-texture'
        
        disp('Computing random tensor field.');
        options.sigma_tensor = getoptions(options, 'lic_regularity', 50*n/256);
        T = compute_tensor_field_random(n,options);
        Flow = perform_tensor_decomp(T); % extract eigenfield.
        options.isoriented = 0; % no orientation in streamlines
        % initial texture
        lic_width = getoptions(options, 'lic_width', 0);
        M0 = perform_blurring(randn(n),lic_width);
        M0 = perform_histogram_equalization( M0, 'linear');
        options.histogram = 'linear';
        options.dt = 0.4;
        options.M0 = M0;
        options.verb = 1;
        options.flow_correction = 1;
        options.niter_lic = 3;
        w = 30;
        M = perform_lic(Flow, w, options);
               
    case 'square_texture'
        M = load_image('square',n);
        M = rescale(M);
        % make a texture patch
        x = linspace(0,1,n);
        [Y,X] = meshgrid(x,x);
        theta = pi/3;
        x = cos(theta)*X + sin(theta)*Y;
        c = [0.3,0.4]; r = 0.2;
        I = find( (X-c(1)).^2 + (Y-c(2)).^2 < r^2 );
        eta = 3/n; lambda = 0.3;
        M(I) = M(I) + lambda * sin( x(I) * 2*pi / eta ); 
        
    case 'tv-image'
        M = rand(n);
        tau = compute_total_variation(M);        
        options.niter = 400;
        [M,err_tv,err_l2] = perform_tv_projection(M,tau/1000,options);
        M = perform_histogram_equalization(M,'linear');
        
    case 'oscillatory_texture'
        x = linspace(0,1,n);
        [Y,X] = meshgrid(x,x);
        theta = pi/3;
        x = cos(theta)*X + sin(theta)*Y;
        c = [0.3,0.4]; r = 0.2;
        I = find( (X-c(1)).^2 + (Y-c(2)).^2 < r^2 );
        eta = 3/n; lambda = 0.3;
        M = sin( x * 2*pi / eta ); 
        
    case {'line', 'line_vertical', 'line_horizontal', 'line_diagonal'}
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        if gamma~=Inf
            M = (X-eta) - gamma*Y < 0;
        else
            M = (Y-eta) < 0;
        end


    case 'line-windowed'
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        eta = .3; 
        gamma = getoptions(options, 'gamma', pi/10);
        parabola = getoptions(options, 'parabola', 0);
        M = (X-eta) - gamma*Y - parabola*Y.^2 < 0;
        f = sin( pi*x ).^2;
        M = M .* ( f'*f );
        
    case 'grating'
        x = linspace(-1,1,n);
        [Y,X] = meshgrid(x,x);
        theta = getoptions(options, 'theta', .2);
        freq = getoptions(options, 'freq', .2);
        X = cos(theta)*X + sin(theta)*Y;
        M = sin(2*pi*X/freq);
        
    case 'disk'
        if ~isfield( options, 'radius' )
            radius = 0.35;
        end
        if ~isfield( options, 'center' )
            center = [0.5, 0.5];    % center of the circle
        end
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        M = (X-center(1)).^2 + (Y-center(2)).^2 < radius^2;
        
        
    case 'diskregular'
        M = rescale(load_image('disk',n,options));
        if not(isfield(options, 'alpha'))
            options.alpha = 3;
        end
        S = load_image('fnoise',n,options);
        M = M + rescale(S,-0.3,0.3); 
        
    case 'quarterdisk'
        if ~isfield( options, 'radius' )
            radius = 0.95;
        end
        if ~isfield( options, 'center' )
            center = -[0.1, 0.1];    % center of the circle
        end
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        M = (X-center(1)).^2 + (Y-center(2)).^2 < radius^2;
        
    case 'fading_contour'
        if ~isfield( options, 'radius' )
            radius = 0.95;
        end
        if ~isfield( options, 'center' )
            center = -[0.1, 0.1];    % center of the circle
        end
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        M = (X-center(1)).^2 + (Y-center(2)).^2 < radius^2;
        theta = 2/pi*atan2(Y,X);
        h = 0.5;
        M = exp(-(1-theta).^2/h^2).*M;
        
    case '3contours'        
        radius = 1.3;
        center = [-1, 1];
        radius1 = 0.8;
        center1 = [0, 0];
        x = 0:1/(n-1):1;
        [Y,X] = meshgrid(x,x);
        f1 = (X-center(1)).^2 + (Y-center(2)).^2 < radius^2;
        f2 = (X-center1(1)).^2 + (Y-center1(2)).^2 < radius1^2;
        M = f1 + 0.5*f2.*(1-f1);
        
    case 'line_circle'

        gamma = 1/sqrt(2);

        x = linspace(-1,1,n);
        [Y,X] = meshgrid(x,x);
        M1 = double( X>gamma*Y+0.25 );
        M2 = X.^2 + Y.^2 < 0.6^2;
        M = 20 + max(0.5*M1,M2) * 216;
        
    case 'fnoise'
        
        % generate an image M whose Fourier spectrum amplitude is 
        %   |M^(omega)| = 1/f^{omega}
        alpha = getoptions(options, 'alpha', 1);
        M = gen_noisy_image(n,alpha);
        
        
    case 'gaussiannoise'
        % generate an image of filtered noise with gaussian
        sigma = getoptions(options, 'sigma', 10);
        M = randn(n);
        m = 51;
        h = compute_gaussian_filter([m m],sigma/(4*n),[n n]);
        M = perform_convolution(M,h);
        return;
    
    case {'bwhorizontal','bwvertical','bwcircle'}
        
        [Y,X] = meshgrid(0:n-1,0:n-1);
        if strcmp(type, 'bwhorizontal')
            d = X;
        elseif strcmp(type, 'bwvertical')
            d = Y;
        elseif strcmp(type, 'bwcircle')
            d = sqrt( (X-(n-1)/2).^2 + (Y-(n-1)/2).^2 );
        end
        if isfield(options, 'stripe_width')
            stripe_width = options.stripe_width;
        else
            stripe_width = 5;
        end
        if isfield(options, 'black_prop')
            black_prop = options.black_prop;
        else
            black_prop = 0.5;
        end
        M = double( mod( d/(2*stripe_width),1 )>=black_prop );
        
    case 'parabola'
        
        % curvature
        c = getoptions(c, 'c', .1);
        % angle
        theta = getoptions(options, 'theta',  pi/sqrt(2));
        x = -0.5:1/(n-1):0.5;
        [Y,X] = meshgrid(x,x);
        Xs = X*cos(theta) + Y*sin(theta);
        Y =-X*sin(theta) + Y*cos(theta); X = Xs;
        M = Y>c*X.^2; 
        
    case 'sin'
        
        [Y,X] = meshgrid(-1:2/(n-1):1, -1:2/(n-1):1);
        M = Y >= 0.6*cos(pi*X);
        M = double(M);
        
    case 'circ_oscil'

        x = linspace(-1,1,n);
        [Y,X] = meshgrid(x,x);
        R = sqrt(X.^2+Y.^2);
        M = cos(R.^3*50);

    case 'phantom'
        
        M = phantom(n);
        
    case 'periodic_bumps'
        nbr_periods = getoptions(options, 'nbr_periods', 8);
        theta = getoptions(options, 'theta', 1/sqrt(2));
        skew = getoptions(options, 'skew', 1/sqrt(2) );        
        A = [cos(theta), -sin(theta); sin(theta), cos(theta)];
        B = [1 skew; 0 1];
        T = B*A;
        x = (0:n-1)*2*pi*nbr_periods/(n-1);
        [Y,X] = meshgrid(x,x);
        pos = [X(:)'; Y(:)'];
        pos = T*pos;
        X = reshape(pos(1,:), n,n);
        Y = reshape(pos(2,:), n,n);
        M = cos(X).*sin(Y);      
        
    case 'noise'
        sigma = getoptions(options, 'sigma', 1);
        M = randn(n) * sigma;
        
    case 'disk-corner'
        x = linspace(0,1,n);
        [Y,X] = meshgrid(x,x);
        rho = .3; eta = .1;
        M1 = rho*X+eta<Y;
        c = [0 .2]; r = .85;
        d = (X-c(1)).^2 + (Y-c(2)).^2;
        M2 = d<r^2;
        M = M1.*M2;
        
    otherwise
        ext = {'gif', 'png', 'jpg', 'bmp', 'tiff', 'pgm', 'ppm'};
        for i=1:length(ext)
            name = [type '.' ext{i}];
            if( exist(name) )
                M = imread( name );
                M = double(M);
                if not(isempty(n)) && (n~=size(M, 1) || n~=size(M, 2)) && nargin>=2
                    M = image_resize(M,n,n);
                end
                if strcmp(type, 'peppers-bw')
                    M(:,1) = M(:,2);
                    M(1,:) = M(2,:);
                end
                if sigma>0
                    M = perform_blurring(M,sigma);
                end
                return;
            end
        end
        error( ['Image ' type ' does not exists.'] );
end

M = double(M);

if sigma>0
    M = perform_blurring(M,sigma);
end

M = rescale(M);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = create_letter(a, r, n)

c = 0.2;
p1 = [c;c];
p2 = [c; 1-c];
p3 = [1-c; 1-c];
p4 = [1-c; c];
p4 = [1-c; c];
pc = [0.5;0.5];
pu = [0.5; c];

switch a
    case 'x'
        point_list = { [p1 p3] [p2 p4] };
    case 'z'
        point_list = { [p2 p3 p1 p4] };
    case 'v'
        point_list = { [p2 pu p3] };
    case 'y'
        point_list = { [p2 pc pu] [pc p3] };
        
        
end
% fit image
for i=1:length(point_list)
    a = point_list{i}(2:-1:1,:);
    a(1,:) = 1-a(1,:);
    point_list{i} = round( a*(n-1)+1 );
end
M = draw_polygons(zeros(n),r,point_list);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sk = draw_polygons(mask,r,point_list)

sk = mask*0;
for i=1:length(point_list)
    pl = point_list{i};
    for k=2:length(pl)
        sk = draw_line(sk,pl(1,k-1),pl(2,k-1),pl(1,k),pl(2,k),r);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sk = draw_line(sk,x1,y1,x2,y2,r)


n = size(sk,1);
[Y,X] = meshgrid(1:n,1:n);
q = 100;
t = linspace(0,1,q);
x = x1*t+x2*(1-t); y = y1*t+y2*(1-t);
if r==0
    x = round( x ); y = round( y );
    sk( x+(y-1)*n ) = 1;
else
    for k=1:q
        I = find((X-x(k)).^2 + (Y-y(k)).^2 <= r^2 );
        sk(I) = 1;
    end
end



function M = gen_noisy_image(n,alpha)

% gen_noisy_image - generate a noisy cloud-like image.
%
%   M = gen_noisy_image(n,alpha);
%
% generate an image M whose Fourier spectrum amplitude is 
%   |M^(omega)| = 1/f^{omega}
%
%   Copyright (c) 2004 Gabriel Peyr?

if nargin<1
    n = 128;
end
if nargin<2
    alpha = 1.5;
end

if mod(n(1),2)==0
    x = -n/2:n/2-1;
else
    x = -(n-1)/2:(n-1)/2;
end

[Y,X] = meshgrid(x,x);
d = sqrt(X.^2 + Y.^2) + 0.1;
f = rand(n)*2*pi;

M = (d.^(-alpha)) .* exp(f*1i);
% M = real(ifft2(fftshift(M)));

M = ifftshift(M);
M = real( ifft2(M) );


function y = gen_signal_2d(n,alpha)

% gen_signal_2d -  generate a 2D C^\alpha signal of length n x n.
%   gen_signal_2d(n,alpha) generate a 2D signal C^alpha. 
%
%   The signal is scale in [0,1].
%   
%   Copyright (c) 2003 Gabriel Peyr?



% new new method

[Y,X] = meshgrid(0:n-1, 0:n-1);

A = X+Y+1;
B = X-Y+n+1;

a = gen_signal(2*n+1, alpha);
b = gen_signal(2*n+1, alpha);
y = a(A).*b(B);
% M = a(1:n)*b(1:n)';

return;


% new method
h = (-n/2+1):(n/2); h(n/2)=1;
[X,Y] = meshgrid(h,h);
h = sqrt(X.^2+Y.^2+1).^(-alpha-1/2);
h = h .* exp( 2i*pi*rand(n,n) );
h = fftshift(h);
y = real( ifft2(h) );

m1 = min(min(y));
m2 = max(max(y));
y = (y-m1)/(m2-m1);

return;

%% old code

y = rand(n,n); 
y = y - mean(mean(y));
for i=1:alpha
    y = cumsum(cumsum(y)')';
    y = y - mean(mean(y));
end
m1 = min(min(y));
m2 = max(max(y));
y = (y-m1)/(m2-m1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = draw_rectangle(r,n)

x = linspace(0,1,n);
[Y,X] = meshgrid(x,x);
M = double( (X>=r(1)) & (X<=r(3)) & (Y>=r(2)) & (Y<=r(4)) ) ;
