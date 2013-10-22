function x = perform_wavelet_transf(x, Jmin, dir, options)

% perform_wavelet_transf - peform fast lifting transform
%
%   y = perform_wavelet_transf(x, Jmin, dir, options);
%
%   Implement 1D and 2D symmetric wavelets with symmetric boundary treatements, using
%   a lifting implementation.
%
%   h = options.filter gives the coefficients of the lifting filter.
%   You can use h='linear' or h='7-9' to select automatically biorthogonal
%   transform with 2 and 4 vanishing moments.
%
%   You can set options.ti=1 to compute a translation invariant wavelet
%   transform.
%
%   You can set options.separable=1 to compute a separable 2D wavelet
%   transform.
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
h = getoptions(options, 'filter', '9-7');
separable = getoptions(options, 'separable', 0);


% detect dimensionality 
% d = getoptions(options,'ndims',-1);
d = -1;
if isfield(options, 'ndims')
    d = options.ndims;
end
if d<0
    if size(x,1)==1 || size(x,2)==1
        d = 1;
    else
        d = 2;
    end
end


if isstr(h)
	% P/U/P/U/etc the last coefficient is scaling
   switch lower(h)
       case 'haar'
           if d==1 || separable==0
               x = perform_haar_transf(x, Jmin, dir, options);
               return;
           end
       case {'linear' '5-3'}
           h = [1/2 1/4 sqrt(2)];
       case {'9-7' '7-9'}
           h = [1.586134342 -.05298011854 -.8829110762 .4435068522 1.149604398];
       otherwise 
           error('Unknown filter');
   end
end


if d==1 && size(x,1)==1
    warning('For 1D transform, the vector should be a column vector.');
    x = x(:);
end

ti = getoptions(options, 'ti', 0);

if d==2 && separable==1
    options.ti = 0;
    if ti==1
        warning('Separable does not works for translation invariant transform');
    end
    % perform a separable wavelet transform
    n = size(x,1);
    if dir==1
        for i=1:n
            x(:,i) = perform_wavelet_transf(x(:,i),Jmin,dir,options);
        end
        for i=1:n
            x(i,:) = perform_wavelet_transf(x(i,:)',Jmin,dir,options)';
        end
    else
        for i=1:n
            x(i,:) = perform_wavelet_transf(x(i,:)',Jmin,dir,options)';
        end
        for i=1:n
            x(:,i) = perform_wavelet_transf(x(:,i),Jmin,dir,options);
        end
    end
    return;
end


% number of lifting steps
n = size(x,1);
m = (length(h)-1)/2;
Jmax = log2(n)-1;
jlist = Jmax:-1:Jmin;
if dir==-1
    jlist = jlist(end:-1:1);
end

if ti==0
    % subsampled
    for j=jlist
        if d==1
            x(1:2^(j+1),:) = lifting_step( x(1:2^(j+1),:), h, dir);
        else
            x(1:2^(j+1),1:2^(j+1)) = lifting_step( x(1:2^(j+1),1:2^(j+1)), h, dir);
            x(1:2^(j+1),1:2^(j+1)) = lifting_step( x(1:2^(j+1),1:2^(j+1))', h, dir)';
        end
    end
else
    % TI
    nJ = Jmax-Jmin+1;
    if dir==1 && d==1
        x = repmat(x,[1 1 nJ+1]);
    elseif dir==1 && d==2
        x = repmat(x,[1 1 3*nJ+1]);
    end
    for j=jlist
        dist = 2^(Jmax-j);
        if d==1
            if dir==1
                x(:,:,[1 j-Jmin+2]) = lifting_step_ti( x(:,:,1), h, dir, dist);
            else
                x(:,:,1) = lifting_step_ti( x(:,:,[1 j-Jmin+2]), h, dir, dist);                
            end
        else
            dj = 3*(j-Jmin);
            if dir==1
                x(:,:,[1 dj+2]) = lifting_step_ti( x(:,:,1), h, dir, dist);
                
                x(:,:,[1 dj+3]) = lifting_step_ti( x(:,:,1)', h, dir, dist);
                x(:,:,1) = x(:,:,1)'; x(:,:,dj+3) = x(:,:,dj+3)';
                
                x(:,:,[dj+[2 4]]) = lifting_step_ti( x(:,:,dj+2)', h, dir, dist);
                x(:,:,dj+2) = x(:,:,dj+2)';
                x(:,:,dj+4) = x(:,:,dj+4)';                
            else                
                x(:,:,dj+2) = x(:,:,dj+2)';
                x(:,:,dj+4) = x(:,:,dj+4)';
                x(:,:,dj+2) = lifting_step_ti( x(:,:,[dj+[2 4]]), h, dir, dist)';

                x(:,:,1) = x(:,:,1)'; 
                x(:,:,dj+3) = x(:,:,dj+3)';
                x(:,:,1) = lifting_step_ti( x(:,:,[1 dj+3]), h, dir, dist)';
                
                x(:,:,1) = lifting_step_ti( x(:,:,[1 dj+2]), h, dir, dist);
           end            
        end
    end
    if dir==-1
        x = x(:,:,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = lifting_step(x, h, dir)

% number of lifting steps
m = (length(h)-1)/2;

if dir==1
    % split
    d = x(2:2:end,:);
    x = x(1:2:end,:);
    for i=1:m
        d = d - h(2*i-1) * ( x + x([2:end end],:) );
        x = x + h(2*i  ) * ( d + d([1 1:end-1],:) );
    end
    x = [x*h(end);d/h(end)];
else
    % retrieve detail coefs
    d = x(end/2+1:end,:)*h(end);
    x = x(1:end/2,:)/h(end);
    for i=m:-1:1
        x = x - h(2*i  ) * ( d + d([1 1:end-1],:) );
        d = d + h(2*i-1) * ( x + x([2:end end],:) );
    end
    % merge
    x1 = zeros(size(x,1)*2,size(x,2));
    x1(1:2:end,:) = x;
    x1(2:2:end,:) = d;
    x = x1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = lifting_step_ti(x, h, dir, dist)

% number of lifting steps
m = (length(h)-1)/2;
n = size(x,1);

s1 = (1:n)+dist;
s2 = (1:n)-dist;
% boundary conditions
s1(s1>n) = 2*n-s1(s1>n); s1(s1<1) = 2-s1(s1<1);
s2(s2>n) = 2*n-s2(s2>n); s2(s2<1) = 2-s2(s2<1);

% s1 = [2 1:n-1]; s2 = [2:n n-1];

if dir==1
    % split
    d = x;
    for i=1:m
        d = d - h(2*i-1) * ( x(s1,:) + x(s2,:) );
        x = x + h(2*i  ) * ( d(s1,:) + d(s2,:) );
    end
    x = cat(3,x*h(end),d/h(end));
else
    % retrieve detail coefs
    d = x(:,:,2)*h(end);
    x = x(:,:,1)/h(end);
    for i=m:-1:1
        x = x - h(2*i  ) * ( d(s1,:) + d(s2,:) );
        d = d + h(2*i-1) * ( x(s1,:) + x(s2,:) );
    end
    % merge
    x = (x+d)/2;
end