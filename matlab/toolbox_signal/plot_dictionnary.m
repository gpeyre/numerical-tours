function H = plot_dictionnary(D,X,nb, options)

% plot_dictionnary - display a dictionary of images
%
%   H = plot_dictionnary(D,X,nb, options);
%
%   options.ndim gives dimensionality (1D/2D)
%   options.s gives number of channel
%
%   D is a (n,p) dictionary of p atoms of dimension n.
%
%   each atom can be a 2D image (ndim==2) of size (s,sqrt(n),sqrt(n))
%   or a 1D signal (ndim==1) of size n
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<2
    X = [];
end
if nargin<3
    nb = 10;
end
if length(nb)==1
    nb = [nb nb];
end


options.null = 0;
ndim = getoptions(options, 'ndim', 2);
s = getoptions(options, 's', 1);

n = size(D,1);
K = size(D,2);

w1 = sqrt(n/s);
pad = 2;

col = [1 1 1];

if ~isempty(X)
    d = sum( abs(X') );
    [tmp,I] = sort(d); I=I(end:-1:1);
else
%    I = randperm(K);
    I = round( linspace(1,K,prod(nb)) );
end

if ndim==1
    % display for 1D dictionary 
    clf;
    k = 0;
    hold on;
    for i=1:nb(1)
        for j=1:nb(2)
            k = k+1;
            if k<=length(I)
                subplot(nb(1), nb(2), k);
                plot(  D(:,I(k)) ); axis tight;
            end
        end
    end
    hold off;
    return;
end

% nb = min(15, floor(sqrt(K)));
H = repmat( reshape(col,[1 1 3]), [nb*(w1+pad) 1] );

vmax = max( max( abs( D(:,I(1:prod(nb))) ) ) );

normalization = getoptions(options, 'normalization', 'rescale');

k = 0;
for i=1:nb(1)
    for j=1:nb(2)
        k = k+1;
        if k<=length(I) % rescale
            v = D(:,I(k));
            if strcmp(normalization, 'clamp')
                v = clamp( 3*v/vmax, -1,1 );
                v = (v+1)/2;
            else
                v = rescale(v);
            end
            v = reshape(v,w1,w1,s);
            selx = (i-1)*(w1+pad)+1:(i-1)*(w1+pad)+w1;
            sely = (j-1)*(w1+pad)+1:(j-1)*(w1+pad)+w1;
            if s==1
                v = repmat(v,[1 1 3]);
            end
            H(selx,sely,:) = v;
        end
    end
end

H(end-pad+1:end,:,:) = [];
H(:,end-pad+1:end,:) = [];

if nargout==0
    clf;
    imagesc(H); axis image; axis off;
end