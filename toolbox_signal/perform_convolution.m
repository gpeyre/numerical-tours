function y = perform_convolution(x,h, bound)

% perform_convolution - compute convolution with centered filter.
%
%   y = perform_convolution(x,h,bound);
%
%   The filter 'h' is centred at 0 for odd
%   length of the filter, and at 1/2 otherwise.
%
%   This works either for 1D or 2D convolution.
%   For 2D the matrix have to be square.
%
%   'bound' is either 'per' (periodic extension) 
%   or 'sym' (symmetric extension).
%
%   Copyright (c) 2004 Gabriel Peyre

if nargin<3
    bound = 'sym';
end

if isstruct(bound)
    bound = getoptions(bound, 'bound', 'sym');
end

if not(strcmp(bound, 'sym')) && not(strcmp(bound, 'per'))
    error('bound should be sym or per');
end

if iscell(x)
    for i=1:length(x)
        y{i} = perform_convolution(x{i},h, bound);
    end
    return;
end

if size(x,3)>1 && size(x,3)<4
    % for color images
    y = x;
    for i=1:size(x,3)
        y(:,:,i) = perform_convolution(x(:,:,i),h, bound);
    end
    return;
end

if size(x,3)>1
    error('Not yet implemented for 3D array, use smooth3 instead.');
end

n = size(x);
p = size(h);

bound = lower(bound);

nd = ndims(x);
if size(x,1)==1 || size(x,2)==1
    nd = 1;
end
if nd==1 
    n = length(x);
    p = length(h);
end

if strcmp(bound, 'sym')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % symmetric boundary conditions
    
    d1 = floor( p/2 );  % padding before
    d2 = p-d1-1;            % padding after
        
    if nd==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% 1D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = x(:); h = h(:);
        xx = [ x(d1:-1:1); x; x(end:-1:end-d2+1) ];
        y = conv(xx,h);
        y = y(p:end-p+1);
    elseif nd==2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % double symmetry
        xx = x;
        xx = [ xx(d1(1):-1:1,:); xx; xx(end:-1:end-d2(1)+1,:) ];
        xx = [ xx(:,d1(2):-1:1), xx, xx(:,end:-1:end-d2(2)+1) ];
        
        y = conv2(xx,h);
        y = y( (2*d1(1)+1):(2*d1(1)+n(1)), (2*d1(2)+1):(2*d1(2)+n(2)) );
    end

else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % periodic boundary conditions
    
    if p>n
        error('h filter should be shorter than x.');
    end
    d = floor((p-1)/2);
    if nd==1    
        x = x(:); h = h(:);
        h = [ h(d+1:end); zeros(n-p,1); h(1:d) ];
        y = real( ifft( fft(x).*fft(h) ) );
    else
        h = [ h(d(1)+1:end,:); zeros( n(1)-p(1),p(2) ); h( 1:d(1),: ) ];
        h = [ h(:,d(2)+1:end), zeros(n(1),n(2)-p(2)), h(:,1:d(2)) ];
        y = real( ifft2( fft2(x).*fft2(h) ) );
    end
    
end