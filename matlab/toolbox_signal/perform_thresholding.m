function y = perform_thresholding(x, t, type, options)

% perform_thresholding - perform hard or soft thresholding
%
%   y = perform_thresholding(x, t, type, options);
%
%   t is the threshold
%   type is either 'hard' or 'soft' or 'semisoft' or 'strict' or 'block'.
%
%   works also for complex data, and for cell arrays.
%
%   if type is 'strict' then it keeps the t largest entry in each
%   column of x.
%
%   for block thresholding, you can set options.block_size to determing the
%   block size.
%
%   Copyright (c) 2006 Gabriel Peyre

options.null = 0;
if nargin==3 && isstruct(type)
    options = type; 
    type = getoptions(options, 'type', 'hard');
end
if nargin<3
    type = 'hard';
end

if iscell(x) && not(strcmp(type, 'freeze'))
    % for cell arrays
    for i=1:size(x,1)
        for j=1:size(x,2)
            y{i,j} = perform_thresholding(x{i,j},t, type);
        end
    end
    return;
end

switch lower(type)
    case {'hard', ''}
        y = perform_hard_thresholding(x,t);
    case 'soft'
        y = perform_soft_thresholding(x,t);
    case 'semisoft'
        y = perform_semisoft_thresholding(x,t);
    case 'strict'
        y = perform_strict_thresholding(x,t);
    case 'largest'
        y = perform_largest_thresh(x,t);
    case {'block' 'block-hard' 'block-soft'}
        if strcmp(type, 'block')
            type = 'block-soft';
        end
        bs = getoptions(options, 'block_size', 4);
        y = perform_block_thresholding(x,t,bs, type);
    case 'quantize'
        y = perform_quant_thresholding(x,t);
    case 'freeze'
        y = perform_freeze_thresholding(x,t);
    case 'soft-multichannel'
        y = perform_softm_thresholding(x,t);
    case {'stein' 'soft-quad'}
        y = perform_stein_thresh(x,t);
    otherwise
        error('Unkwnown thresholding type.');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = perform_strict_thresholding(X,s)
%% keep only the s largest coefficients in each column of X

v = sort(abs(X)); v = v(end:-1:1,:);
v = v(round(s),:);
v = repmat(v, [size(X,1) 1]);
X = X .* (abs(X)>=v);

function X = perform_largest_thresh(X,s)

v = sort(abs(X(:))); v = v(end:-1:1,:);
v = v(round(s));
X = X .* (abs(X)>=v);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_hard_thresholding(x,t)

t = t(1);
y = x .* (abs(x) > t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_quant_thresholding(x,t)

t = t(1);
y = floor(abs(x/t)).*sign(x);
y = sign(y) .* (abs(y)+.5) * t;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_stein_thresh(x,t)

t = t(1);
y = x .* max( 1-t^2 ./ abs(x).^2, 0 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_softm_thresholding(x,t)

% dimension of thresholding
dm = nb_dims(x);

% soft multichannel
d = sqrt( sum( abs(x).^2, dm) );
d = repmat(d, [ones(dm-1,1); size(x,dm)]);
y = max( 1-t./d, 0 ) .* x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_soft_thresholding(x,t)

if not(isreal(x))
    % complex threshold
    d = abs(x);
    d(d<eps) = 1;
    y = x./d .* perform_soft_thresholding(d,t);
    return;
end

%if nb_dims(x)==nb_dims(t) && size(t)~=size(x)
    t = t(1);
%end
s = abs(x) - t;
s = (s + abs(s))/2;
y = sign(x).*s;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_freeze_thresholding(x,t)

if iscell(x)
    if length(x)~=length(t)
        error('Vector should have same length');
    end
    for i=1:length(x)
        y{i} = perform_freeze_thresholding(x{i},t{i});
    end
    return;
end

y = x;
y(t==0) = 0;
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_semisoft_thresholding(x,t)

if length(t)==1
    t = [t 2*t];
end
t = sort(t);

y = x;
y(abs(x)<t(1)) = 0;
I = find(abs(x)>=t(1) & abs(x)<t(2));
y( I ) = sign(x(I)) .* t(2)/(t(2)-t(1)) .* (abs(x(I))-t(1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = perform_block_thresholding(x,t,bs, type)

n = size(x,1);
if nb_dims(x)==2
    p = size(x,2);
    if length(bs)==1
        bs = [bs bs];
    end
    % compute indexing
    [dX,dY,X,Y] = ndgrid(0:bs(1)-1,0:bs(2)-1,1:bs(1):n-bs(1)+1,1:bs(2):p-bs(2)+1);
    I = X+dX + (Y+dY-1)*n;
    % reshape as block
    H = x(I);
    % threshold
    v = mean(mean(abs(H).^2,1),2); v = max(v,1e-15);
    if strcmp(type, 'block-soft')
        H = repmat(max(1-t^2./v,0),[bs(1) bs(2) 1 1]) .* H;
    else
        H = repmat(v>t^2,[bs bs 1 1]) .* H;
    end
elseif nb_dims(x)==1
    %% 1D %%
    [dX,X] = meshgrid(0:bs-1,1:bs:n-bs+1);
    I = X+dX;
    % reshape as block
    H = x(I);
    % threshold
    v = mean(H.^2); v = max(v,1e-15);
    if strcmp(type, 'block-soft')
        H = repmat(max(1-t^2./v,0),[bs 1]) .* H;
    else
        H = double( repmat(v>t^2,[bs 1]) ) .* H;
    end
elseif nb_dims(x)==3
    for i=1:size(x,3)
        y(:,:,i) = perform_block_thresholding(x(:,:,i),t,bs, type);
    end
    return;
else
    error('Wrong size');
end
% reconstruct
y = x; y(I) = H;