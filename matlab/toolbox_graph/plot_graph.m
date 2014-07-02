function h = plot_graph(A,xy, options)

% plot_graph - display a 2D or 3D graph.
%
% plot_graph(A,xy, options);
%
%   options.col set the display (e.g. 'k.-')
%
%   Copyright (c) 2006 Gabriel Peyre

if size(xy,1)>size(xy,2)
    xy = xy';
end

if nargin<3
    options.null = 0;
end

if not(isstruct(options))
    col = options; clear options;
    options.null = 0;
else
    if isfield(options, 'col')
        col = options.col;
    else
        col = 'k.-';
    end
end
lw = getoptions(options, 'lw', 2);
ps = getoptions(options, 'ps', 20);


B = full(A); B = B(:); B = B( B~=0 );
isvarying = std(B)>0;

if size(xy,1)==2
    if ~isstr(col)
        col = [];
    end
    % 2D display
    h = gplotvarying(A,xy, col);
    % h = gplot(A,xy', col);
    axis tight; axis equal; axis off;
elseif size(xy,1)==3
    % 3D display
    if isstr(col)
        if ~isvarying
            h = gplot3(A,xy', col);
        else
            h = gplotvarying3(A,xy, col);
        end
    else
        hold on;
        if ~isvarying
            h = gplot3(A,xy', col);
        else
            h = gplotvarying3(A,xy, col);
        end
        plot_scattered(xy, col);
        hold off;
        view(3);
    end
    axis off;
    cameramenu;
else
    error('Works only for 2D and 3D graphs');
end

set(h, 'LineWidth', lw);
set(h, 'MarkerSize', ps);

function h = gplotvarying(A,xy,col)


[i,j,s] = find(sparse(A));
I = find(i<=j); 
i = i(I); j = j(I); s = s(I);

x = [xy(1,i);xy(1,j)];
y = [xy(2,i);xy(2,j)];
h = plot(x,y,col);


function h = gplotvarying3(A,xy,col)

[i,j,s] = find(sparse(A));
I = find(i<=j); 
i = i(I); j = j(I); s = s(I);

x = [xy(1,i);xy(1,j)];
y = [xy(2,i);xy(2,j)];
y = [xy(3,i);xy(3,j)];
h = plot3(x,y,z,col);

function h = gplot3(A,xy,lc)


[i,j] = find(A);
[ignore, p] = sort(max(i,j));
i = i(p);
j = j(p);

% Create a long, NaN-separated list of line segments,
% rather than individual segments.

X = [ xy(i,1) xy(j,1) repmat(NaN,size(i))]';
Y = [ xy(i,2) xy(j,2) repmat(NaN,size(i))]';
Z = [ xy(i,3) xy(j,3) repmat(NaN,size(i))]';
X = X(:);
Y = Y(:);
Z = Z(:);

if nargin<3,
    h = plot3(X, Y, Z)
else
    h = plot3(X, Y, Z, lc);
end