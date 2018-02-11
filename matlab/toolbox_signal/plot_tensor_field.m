function h = plot_tensor_field(H, M, options)

% plot_tensor_field - display a tensor field
%
%   h = plot_tensor_field(H, M, options);
%
%   options.sub controls sub-sampling
%   options.color controls color
%
%   Copyright (c) 2006 Gabriel Peyre

if nargin<3
    options.null = 0;
end
if not( isstruct(options) )
    sub = options;
    clear options;
    options.sub = sub;
end

% sub = getoptions(options, 'sub', 1);
sub = getoptions(options, 'sub', round(size(H,1)/30) );
color = getoptions(options, 'color', 'r');

if nargin<2
    M = [];
end

if not(isempty(M)) && size(M,3)==1
    M = repmat(M, [1 1 3]); % ensure B&W image
end

if size(H,3)==3 && size(H,4)==1
    H = cat(3, H(:,:,1), H(:,:,3), H(:,:,3), H(:,:,2) );
    H = reshape(H, size(H,1), size(H,2), 2, 2);
    if 0
    % flip the main eigen-axes
    [e1,e2,l1,l2] = perform_tensor_decomp(H);
    H = perform_tensor_recomp(e2,e1,l1,l2);
    end
    h = plot_tensor_field(H, M, sub);
    return;
end

% swap X and Y axis
%%% TODO
a = H(:,:,2,2);
H(:,:,2,2) = H(:,:,1,1);
H(:,:,1,1) = a;

hold on;
if ~isempty(M)
    imagesc(rescale(M)); drawnow;
end
h = fn_tensordisplay(H(:,:,1,1),H(:,:,1,2), H(:,:,2,2), 'sub', sub, 'color', color);
axis image; axis off;
colormap jet(256);
% hold off;


function h = fn_tensordisplay(varargin)

% function h = fn_tensordisplay([X,Y,]Txx,Txy,Tyy[,'sigma',sigma][,'sub',sub][,color][,patch options...]])
% function h = fn_tensordisplay([X,Y,]e[,'sigma',sigma][,'sub',sub][,color][,patch options...]])

% X,Y,Txx,Txy,Tyy
if isstruct(varargin{1}) || isstruct(varargin{3})
    if isstruct(varargin{1}), nextarg=1; else nextarg=3; end
    e = varargin{nextarg};
    Txx = e.ytyt;
    Txy = -e.ytyx;
    Tyy = e.yxyx;
    if nextarg==1
        [nj ni] = size(Txx);
        [X Y] = meshgrid(1:ni,1:nj);
    else
        [X Y] = deal(varargin{1:2});
    end
    nextarg = nextarg+1;
else
    [nj ni] = size(varargin{3});
    if nargin<5 || ischar(varargin{4}) || ischar(varargin{5}) || any(size(varargin{5})~=[nj ni])
        [X Y] = meshgrid(1:ni,1:nj);
        nextarg = 1;
    else
        [X Y] = deal(varargin{1:2});
        nextarg = 3;
    end
    [Txx Txy Tyy] = deal(varargin{nextarg:nextarg+2});
    nextarg = nextarg+3;
end
if any(size(X)==1), [X Y] = meshgrid(X,Y); end
[nj,ni] = size(X);
if any(size(Y)~=[nj ni]) || ...
        any(size(Txx)~=[nj ni]) || any(size(Txy)~=[nj ni]) || any(size(Tyy)~=[nj ni])
    error('Matrices must be same size')
end
% sigma, sub, color
color = 'r';
while nextarg<=nargin
    flag = varargin{nextarg};
    nextarg=nextarg+1;
    if ~ischar(flag), color = flag; continue, end
    switch lower(flag)
        case 'sigma'
            sigma = varargin{nextarg};
            nextarg = nextarg+1;
            switch length(sigma)
                case 1
                    sigmax = sigma;
                    sigmay = sigma;
                case 2
                    sigmax = sigma(1);
                    sigmay = sigma(2);
                otherwise
                    error('sigma definition should entail two values');
            end
            h = fspecial('gaussian',[ceil(2*sigmay) 1],sigmay)*fspecial('gaussian',[1 ceil(2*sigmax)],sigmax);
            Txx = imfilter(Txx,h,'replicate');
            Txy = imfilter(Txy,h,'replicate');
            Tyy = imfilter(Tyy,h,'replicate');
        case 'sub'
            sub = varargin{nextarg};
            nextarg = nextarg+1;
            switch length(sub)
                case 1
                    [x y] = meshgrid(1:sub:ni,1:sub:nj);
                    sub = y+nj*(x-1);
                case 2
                    [x y] = meshgrid(1:sub(1):ni,1:sub(2):nj);
                    sub = y+nj*(x-1);
            end
            X = X(sub); Y = Y(sub);
            Txx = Txx(sub); Txy = Txy(sub); Tyy = Tyy(sub);
            [nj ni] = size(sub);
        case 'color'            
            color = varargin{nextarg};
            nextarg = nextarg+1;
        otherwise
            break
    end
end
% options
options = {varargin{nextarg:end}};


npoints = 50;
theta = (0:npoints-1)*(2*pi/npoints);
circle = [cos(theta) ; sin(theta)];
Tensor = cat(3,Txx,Txy,Txy,Tyy);        % jdisplay x idisplay x tensor
Tensor = reshape(Tensor,2*nj*ni,2);     % (display x 1tensor) x 2tensor
Ellipse = Tensor * circle;              % (display x uv) x npoints
Ellipse = reshape(Ellipse,nj*ni,2,npoints);         % display x uv x npoints
XX = repmat(X(:),1,npoints);                        % display x npoints
YY = repmat(Y(:),1,npoints);                        % display x npoints
U  = reshape(Ellipse(:,1,:),nj*ni,npoints);         % display x npoints
V  = reshape(Ellipse(:,2,:),nj*ni,npoints);         % display x npoints
umax = max(U')'; vmax = max(V')';
umax(umax==0)=1; vmax(vmax==0)=1;
if ni==1, dx=1; else dx = X(1,2)-X(1,1); end
if nj==1, dy=1; else dy = Y(2,1)-Y(1,1); end
fact = min(dx./umax,dy./vmax)*.35;
fact = repmat(fact,1,npoints);
U = XX + fact.*U;
V = YY + fact.*V;

%-----------
MM = mmax(Txx+Tyy);
mm = mmin(Txx+Tyy);
[S1, S2] = size(U);
Colormap = zeros(S2, S1, 3);
Map = colormap(jet(256));
for k =1:npoints
    Colormap(k,:,1) = Map(floor(255*(Txx(:)+Tyy(:)-mm)/(MM-mm)) + 1, 1);
    Colormap(k,:,2) = Map(floor(255*(Txx(:)+Tyy(:)-mm)/(MM-mm)) + 1, 2);
    Colormap(k,:,3) = Map(floor(255*(Txx(:)+Tyy(:)-mm)/(MM-mm)) + 1, 3);
end
%-----------
%h = fill(U',V',color,'EdgeColor',color,options{:});
h = fill(U',V',Colormap,'EdgeColor', 'interp');

axis ij;

if nargout==0, clear h, end


function a=mmax(a)
a = max(a(:));
function a=mmin(a)
a = min(a(:));
